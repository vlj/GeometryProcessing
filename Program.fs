﻿// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

open Emgu.CV
open Emgu.CV.Structure
open System.Drawing
open MathNet.Numerics.LinearAlgebra
open MathNet.Spatial.Euclidean
open System.Numerics
open System.Collections.Generic
open FSharpx.Collections
open MathNet.Numerics



let project [|(p0:Vector3D);(p1:Vector3D);(p2:Vector3D)|] =
    let X = p1 - p0
    let Y = p2 - p0
    let Ex = X.Normalize ()
    let Ez = Y.CrossProduct(Ex).Normalize ()
    let Ey = Ez.CrossProduct(Ex)
    [| Vector2D(p0.DotProduct(Ex), p0.DotProduct(Ey));
       Vector2D(p1.DotProduct(Ex), p1.DotProduct(Ey));
       Vector2D(p2.DotProduct(Ex), p2.DotProduct(Ey))|]

let coeff [|(p1:Vector2D); (p2:Vector2D); (p3:Vector2D)|] =
    let dT = (p1.X * p2.Y - p1.Y * p2.X) + (p2.X * p3.Y - p2.Y * p3.X) + (p3.X * p1.Y - p3.Y * p1.X)
    let s = (sign >> float) dT
    let sdT = Complex(s * (sqrt (s * dT)), 0.)
    let W1 = Complex(p3.X - p2.X, p3.Y - p2.Y) / sdT
    let W2 = Complex(p1.X - p3.X, p1.Y - p3.Y) / sdT
    let W3 = Complex(p2.X - p1.X, p2.Y - p1.Y) / sdT
    (W1, W2, W3)

[<EntryPoint>]
let main argv = 
    use file = CvInvoke.Imread(@"C:\Users\vlj\Desktop\height_map_norway-height-map-aster-30m.png")

    let size = 10
    let vertex_to_idx i j = i + size * j
    let triangle_to_idx i j = i + (size - 1) * j
    let triangle_count = (size - 1) * (size - 1) * 4

    let idx_to_vertex = (function idx -> (idx / size, idx % size))
    let idx_to_triangle = (function idx-> (idx / (size - 1), idx % (size - 1)))
    let get_height x y = let c = file.Bitmap.GetPixel(int(x), int(y)) in float(c.GetBrightness()) * 200.
    let get_points (i,j) = let x, y = float(i) * 100. + 10., float(j) * 100. + 10. in Vector3D(x, y, get_height x y)
    let get_middle (i,j) = let x, y = float(i) * 100. + 60., float(j) * 100. + 60. in Vector3D(x, y, get_height x y)

    let square = Array.init (size * size) (idx_to_vertex >> get_points)
    let middles = Array.init ((size - 1) * (size - 1)) (idx_to_triangle >> get_middle)
    let points = Array.append square middles


    let square_idx_to_triangle idx =
         let (i,j) = idx_to_triangle idx
         let left_up = vertex_to_idx i j
         let left_down = vertex_to_idx (i + 1) j
         let right_down = vertex_to_idx (i + 1) (j + 1)
         let right_up = vertex_to_idx i (j + 1)
         let center = square.Length + (triangle_to_idx i j)
         [|
             [| left_up; left_down; center |];
             [| left_down; right_down; center |];
             [| right_down; right_up; center |];
             [| right_up; left_up; center |]
         |]
    let triangles = Seq.collect square_idx_to_triangle { 0..((size - 1) * (size - 1) - 1) }

    let build_matrix =
        let result = CreateMatrix.Dense<Complex>(triangle_count, points.Length)
        Seq.iteri (fun triangle_idx triangle ->
            let (W1, W2, W3) = coeff ((Array.map (function x -> points.[x]) >> project) triangle)
            let col_idx idx = triangle.[idx]
            result.[triangle_idx, col_idx 0] <- W1
            result.[triangle_idx, col_idx 1] <- W2
            result.[triangle_idx, col_idx 2] <- W3) triangles
        result

    // there's 100 vertex, 10 on top/bottom + 8 on left/right are fixed (36)
    let pinned_vector = CreateVector.Dense<Complex> (4 * (size - 1))
    let start_of_fixed = points.Length - (4 * (size - 1))
    let dic = [|for i in 0..(points.Length - 1) -> i|]
    let next_free_idx = ref 0
    let idx = ref 0
    for i in 0..(size - 1) do
        let bottom = vertex_to_idx i 0
        dic.[bottom] <- (start_of_fixed + !idx)
        pinned_vector.[!idx] <- Complex(0., float(i) * 100.)
        incr idx

    printfn "%A" dic
    for i in 0..(size - 1) do
        let top = vertex_to_idx i (size - 1)
        dic.[top] <- (start_of_fixed + !idx)
        pinned_vector.[!idx] <- Complex(900., float(i) * 100.)
        incr idx

    printfn "%A" dic
    for i in 1..(size - 2) do
        let left = vertex_to_idx 0 i
        dic.[left] <- (start_of_fixed + !idx)
        pinned_vector.[!idx] <- Complex(float(i) * 100., 0.)
        incr idx

    printfn "%A" dic
    for i in 1..(size - 2) do
        let right = vertex_to_idx (size - 1) i
        dic.[right] <- (start_of_fixed + !idx)
        pinned_vector.[!idx] <- Complex(float(i) * 100., 900.)
        incr idx

    printfn "%A" dic
    for i in 1..(size - 2) do
        for j in 1..(size - 2) do
            dic.[vertex_to_idx i j] <- !next_free_idx
            incr next_free_idx
    for i in square.Length..(points.Length - 1) do
        dic.[i] <- !next_free_idx
        incr next_free_idx
    printfn "%A" dic
    build_matrix.PermuteColumns( Permutation(dic))

    let freemat = Seq.fold (fun (s:Matrix<Complex>) i -> s.RemoveColumn i) build_matrix {(points.Length - 1) ..(-1)..start_of_fixed}
    let pinmat = Seq.fold (fun (s:Matrix<Complex>) i -> s.RemoveColumn i) build_matrix {(start_of_fixed - 1)..(-1)..0}

    let bb = -pinmat.Multiply(pinned_vector)
    let uvs = freemat.Solve(bb)

    let to_points (v:Complex) = Point(int(v.Real), int(v.Imaginary))

    let pos2D = [| for i in 0..points.Length - 1 ->
        let new_idx = dic.[i]
        if new_idx < start_of_fixed then
            to_points uvs.[new_idx]
        else
            to_points pinned_vector.[new_idx - start_of_fixed]
        |]

//    Array.iter (function p -> CvInvoke.Circle(file, to_points p, 5, Bgr(0.0, 0.0, 255.0).MCvScalar)) points
    let draw_triangle arr = CvInvoke.Polylines(file, Array.ofSeq arr, true, Bgr(255.0, 0.0, 0.0).MCvScalar)
    let to_pts = Seq.map (Seq.map (function x -> pos2D.[x])) triangles
    Seq.iter draw_triangle to_pts


    CvInvoke.Imshow("", file)
    CvInvoke.WaitKey 0
    0 // return an integer exit code
