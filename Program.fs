// Learn more about F# at http://fsharp.org
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
open MathNet.Numerics.LinearAlgebra.Factorization
open MathNet.Numerics.LinearAlgebra.Matrix


type triangle(indexes : int array) as this =
        let indexes = indexes
        member this.GetVertex (vertexes : Vector3D array) =
            let p0 = vertexes.[indexes.[0]]
            let p1 = vertexes.[indexes.[1]]
            let p2 = vertexes.[indexes.[2]]
            (p0, p1, p2)
        member this.GetLocalBasis (vertexes : Vector3D array) =
            let (p0, p1, p2) = this.GetVertex vertexes
            let X = p1 - p0
            let Y = p2 - p0
            let Ex = X.Normalize ()
            let Ez = Y.CrossProduct(Ex).Normalize ()
            let Ey = Ez.CrossProduct(Ex)
            (Ex, Ey, Ez)
        member this.GetIndexes () = indexes
        member this.GetCotangent (vertexes: Vector3D array) =
            let (A, B, C) = this.GetVertex vertexes
            let AB = B - A
            let BC = C - B
            let CA = A - C
            let a = AB.Length
            let b = BC.Length
            let c = CA.Length
            let args = [| AB, -CA, c * b; -AB, BC, c * a; CA, -BC, a * b |]
            let coss = args |> Array.map (function (e1, e2, denum) -> e1.DotProduct(e2) / denum)
            let sins = args |> Array.map (function (e1, e2, denum) -> e1.CrossProduct(e2).Length / denum)
            Array.map2 (/) coss sins
        member this.project (vertexes: Vector3D array) =
            let (Ex, Ey, Ez) = this.GetLocalBasis vertexes
            let (p0, p1, p2) = this.GetVertex vertexes
            let X = p1 - p0
            let Y = p2 - p0
            [| Vector2D(p0.DotProduct(Ex), p0.DotProduct(Ey));
               Vector2D(p1.DotProduct(Ex), p1.DotProduct(Ey));
               Vector2D(p2.DotProduct(Ex), p2.DotProduct(Ey))|]



let closest_rotation A =
    let v = svd A
    if v.S.[0] * v.S.[1] < 0. then
        (- v.U * v.VT, - v.VT.Transpose () * v.W * v.VT)
    else
        (v.U * v.VT, v.VT.Transpose () * v.W * v.VT)

let mapping (points : Vector3D array) (border_point: IDictionary<int, Vector2D>) =
    let pinned_vector = CreateVector.Dense<Complex> (border_point.Count)
    let start_of_fixed = points.Length - border_point.Count

    let dic = [|for i in 0..(points.Length - 1) -> i|]
    let pinned_idx = ref 0
    for pair in border_point do
        let (idx, v) = pair.Key, pair.Value
        dic.[idx] <- (start_of_fixed + !pinned_idx)
        pinned_vector.[!pinned_idx] <- Complex(v.X, v.Y)
        incr pinned_idx

    let next_free_idx = ref 0
    for i in 0..points.Length - 1 do
        if not (border_point.ContainsKey i) then 
            dic.[i] <- !next_free_idx
            incr next_free_idx
    dic, pinned_vector


let merge_free_and_pinned (free: Vector2D array) (pinned: Vector2D array) (dic: int array) =
    [|
        let total_count = free.Length + pinned.Length
        for i in 0..total_count - 1 do
            let new_idx = dic.[i]
            if new_idx < free.Length then
                yield free.[new_idx]
            else
                yield pinned.[new_idx - free.Length]
    |]


let LSCM (points : Vector3D array) (border_point: IDictionary<int, Vector2D>) (triangles: triangle array) =


    let coeff [|(p1:Vector2D); (p2:Vector2D); (p3:Vector2D)|] =
        let dT = (p1.X * p2.Y - p1.Y * p2.X) + (p2.X * p3.Y - p2.Y * p3.X) + (p3.X * p1.Y - p3.Y * p1.X)
        let s = (sign >> float) dT
        let sdT = Complex(s * (sqrt (s * dT)), 0.)
        let W1 = Complex(p3.X - p2.X, p3.Y - p2.Y) / sdT
        let W2 = Complex(p1.X - p3.X, p1.Y - p3.Y) / sdT
        let W3 = Complex(p2.X - p1.X, p2.Y - p1.Y) / sdT
        (W1, W2, W3)

    let start_of_fixed = points.Length - border_point.Count

    let dic, pinned_vector = mapping points border_point

    let get_freemat_pinmat () =
        let C_matrix = CreateMatrix.Dense<Complex>(triangles.Length, points.Length)
        triangles |> Seq.iteri (fun triangle_idx t ->
            let (W1, W2, W3) = t.project points |> coeff
            let col_idx idx = (t.GetIndexes ()).[idx]
            C_matrix.[triangle_idx, col_idx 0] <- W1
            C_matrix.[triangle_idx, col_idx 1] <- W2
            C_matrix.[triangle_idx, col_idx 2] <- W3)
        C_matrix.PermuteColumns(Permutation(dic))
        let freemat = C_matrix.[0.., ..start_of_fixed - 1]
        let pinmat = C_matrix.[0.., start_of_fixed..]
        freemat, pinmat

    let freemat, pinmat = get_freemat_pinmat()
    let bb = pinmat.Multiply(-pinned_vector)
    let uvs = freemat.Solve(bb)

    let to_points (v:Complex) = Vector2D(v.Real, v.Imaginary)
    let free = [|for v in uvs -> to_points v |]
    let pinned = [|for v in pinned_vector -> to_points v |]
    merge_free_and_pinned free pinned dic

let arap (points : Vector3D array) (border_point: IDictionary<int, Vector2D>) (triangles: triangle array) =
    //let neighbor =
    //    let neighbor_dup = [| for p in points -> ([]:int list) |]
    //    for t in triangles do
    //        let [|i0; i1; i2|] = t.GetIndexes ()
    //        neighbor_dup.[i0] <- i1 :: i2 :: neighbor_dup.[i0]
    //        neighbor_dup.[i1] <- i0 :: i2 :: neighbor_dup.[i1]
    //        neighbor_dup.[i2] <- i1 :: i0 :: neighbor_dup.[i2]
    //    Array.map List.distinct neighbor_dup


    let initial_guess = LSCM points border_point triangles

    let weights =
        let w = [| for p in points -> ([]: (int * Vector<float> * float) list)|]
        let insert_or_sum k v =
            let aux = w.[k] |> List.exists (function (o, _, _) -> let (o2, _, _) = v in o = o2)
            //if aux then failwith "not ok"
            w.[k] <- v :: w.[k]
        for t in triangles do
            let [|i0; i1; i2|] = t.GetIndexes ()
            let cotan = t.GetCotangent points
            let [|x0; x1; x2|] = t.project points
            let x0x1 = (x1 - x0).ToVector ()
            let x1x2 = (x2 - x1).ToVector ()
            let x2x0 = (x0 - x2).ToVector ()
            let zero = CreateVector.Dense<float>(2)
            insert_or_sum i1 (i2, x1x2, cotan.[0] / 2.)
            insert_or_sum i2 (i1, -x1x2, cotan.[0] / 2.)
            insert_or_sum i0 (i2, -x2x0, cotan.[1] / 2.)
            insert_or_sum i2 (i0, x2x0, cotan.[1] / 2.)
            insert_or_sum i0 (i1, x0x1, cotan.[2] / 2.)
            insert_or_sum i1 (i0, -x0x1, cotan.[2] / 2.)

            insert_or_sum i1 (i1, zero, -cotan.[0] / 2.)
            insert_or_sum i2 (i2, zero, -cotan.[0] / 2.)
            insert_or_sum i0 (i0, zero, -cotan.[1] / 2.)
            insert_or_sum i2 (i2, zero, -cotan.[1] / 2.)
            insert_or_sum i0 (i0, zero, -cotan.[2] / 2.)
            insert_or_sum i1 (i1, zero, -cotan.[2] / 2.)
        w

    let energy i (current_uv: Vector2D array) =
        let cross_variance i j (half_edge:Vector<float>) (cot:float) = 
            let uiuj = (current_uv.[i] - current_uv.[j]).ToVector ()
            half_edge.OuterProduct(uiuj).Multiply(cot)
        let add s (j, half_edge, w) = s + (cross_variance i j half_edge w)
        weights.[i] |> List.fold add (CreateMatrix.Dense<float>(2, 2))

    let rotations = [|
        for i in 0..points.Length - 1 ->
            let e = energy i initial_guess
            let (R, _) = closest_rotation e
            R
        |]

    let map, pinned_vector = mapping points border_point
    let start_of_fixed = points.Length - border_point.Count

    let pinned_vector2 = CreateMatrix.Dense(pinned_vector.Count, 2)
    for i in 0..pinned_vector.Count - 1 do
        pinned_vector2.[i, 0] <- pinned_vector.[i].Real
        pinned_vector2.[i, 1] <- pinned_vector.[i].Imaginary

    let WeightMatrix =
        let L = CreateMatrix.Dense<float> (points.Length, points.Length)
        for row in 0..L.RowCount - 1 do
            let assign (col, _, w) =
                L.[row, col] <- L.[row, col] + w
            weights.[row] |> List.iter assign
        L.PermuteColumns(Permutation(map))
        L.PermuteRows(Permutation(map))
        L

    printfn "WeightMatrix : %A" (WeightMatrix.ToString ())

    let laplacian =
        -1. * WeightMatrix.[.. start_of_fixed - 1, .. start_of_fixed - 1]

    printfn "Laplacian : %A" (laplacian.ToString ())

    let rhs =
        let r = WeightMatrix.[.. start_of_fixed - 1, start_of_fixed..] * pinned_vector2
        let add i j (half_edge:Vector<float>) w (R: Matrix<float> array) =
            let M = (R.[i] - R.[j]).TransposeThisAndMultiply(half_edge)
            w / 2. * M
        //for row in 0..points.Length - border_point.Count - 1 do
        //    let i = map.[row]
        //    let rowcontent = weights.[i] |> List.fold (fun s (j, half_edge, w) -> s + (add i j half_edge w rotations)) (CreateVector.Dense<float>(2))
        //    r.SetRow(row, rowcontent)
        r

    let res = laplacian.Solve(rhs)
    let reshaped = [| for row in res.ToRowArrays () -> Vector2D(row.[0], row.[1]) |]
    printfn "reshaped:%A" reshaped
    let pinned = [|for v in pinned_vector -> Vector2D(v.Real, v.Imaginary) |]
    merge_free_and_pinned reshaped pinned map


[<EntryPoint>]
let main argv = 
    use file = CvInvoke.Imread(@"C:\Users\vlj\Desktop\height_map_norway-height-map-aster-30m.png")

    let size = 3
    let vertex_to_idx i j = i + size * j
    let triangle_to_idx i j = i + (size - 1) * j
    let triangle_count = (size - 1) * (size - 1) * 4

    use bmp = file.Bitmap

    let step_x = float(bmp.Width - 10) / float(size)
    let step_y = float(bmp.Height - 10) / float(size)

    let idx_to_vertex = (function idx -> (idx / size, idx % size))
    let idx_to_triangle = (function idx-> (idx / (size - 1), idx % (size - 1)))
    let get_height _x _y =
        let x, y = int(round _x), int(round _y)
        let c = bmp.GetPixel(x, y)
        float(c.GetBrightness()) * 300.
    let get_points (i,j) =
        let x, y = float(i) * step_x, float(j) * step_y
        Vector3D(x, y, get_height x y)
    let get_middle (i,j) =
        let x, y = (float(i) + 0.5) * step_x, (float(j) + 0.5) * step_y
        Vector3D(x, y, get_height x y)

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
             triangle [| left_up; left_down; center |];
             triangle [| left_down; right_down; center |];
             triangle [| right_down; right_up; center |];
             triangle [| right_up; left_up; center |]
         |]
    let triangles = Seq.collect square_idx_to_triangle { 0..((size - 1) * (size - 1) - 1) }

    let pinned_val = dict[
        for i in 0..(size - 1) do
            yield (vertex_to_idx i 0, Vector2D(0., float(i) * step_y))
        for i in 0..(size - 1) do
            yield (vertex_to_idx i (size - 1), Vector2D(step_x * float(size), float(i) * step_y))
        for i in 1..(size - 2) do
            yield (vertex_to_idx 0 i, Vector2D(float(i) * step_x, 0.))
        for i in 1..(size - 2) do
            yield (vertex_to_idx (size - 1) i, Vector2D(float(i) * step_x, step_y * float(size)))
        ]


    let pos2D = arap points pinned_val (Array.ofSeq triangles)

//    Array.iter (function p -> CvInvoke.Circle(file, to_points p, 5, Bgr(0.0, 0.0, 255.0).MCvScalar)) points
    let draw_triangle arr = CvInvoke.Polylines(file, Array.ofSeq arr, true, Bgr(255.0, 0.0, 0.0).MCvScalar)
    let to_pts = Seq.map (function (t:triangle) -> t.GetIndexes () |> Array.map (function x -> let v = pos2D.[x] in Point(int(v.X), int(v.Y)))) triangles
    Seq.iter draw_triangle to_pts


    CvInvoke.Imshow("", file)
    CvInvoke.WaitKey 0
    0 // return an integer exit code
