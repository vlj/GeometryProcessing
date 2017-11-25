module slim

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Factorization
open MathNet.Numerics.LinearAlgebra.Matrix
open triangle
open MathNet.Spatial.Euclidean
open System.Numerics
open System.Collections.Generic
open MathNet.Numerics
open arap

let SLIM (points : Vector3D array) (border_point: IDictionary<int, Vector2D>) (triangles: triangle array) =
    let compute_surface_gradient_matrix (t:triangle) =
        let [|i1;i2;i3|] = t.GetIndexes ()
        let perpx1x3, perpx2x1 = t.getGradient points
        let (Ex, Ey, _) = t.GetLocalBasis points
        let D1 = CreateVector.Dense<float> (points.Length)
        D1.[i2] <- Ex.X * perpx1x3.X + Ex.Y * perpx1x3.Y + Ex.Z * perpx1x3.Z
        D1.[i3] <- Ex.X * perpx2x1.X + Ex.Y * perpx2x1.Y + Ex.Z * perpx2x1.Z
        D1.[i1] <- - Ex.X * (perpx1x3.X - perpx2x1.X) - Ex.Y * (perpx1x3.Y - perpx2x1.Y) - Ex.Z * (perpx1x3.Z - perpx2x1.Z)
        let D2 = CreateVector.Dense<float> (points.Length)
        D2.[i2] <- Ey.X * perpx1x3.X + Ey.Y * perpx1x3.Y + Ey.Z * perpx1x3.Z
        D2.[i3] <- Ey.X * perpx2x1.X + Ey.Y * perpx2x1.Y + Ey.Z * perpx2x1.Z
        D2.[i1] <- - Ey.X * (perpx1x3.X - perpx2x1.X) - Ey.Y * (perpx1x3.Y - perpx2x1.Y) - Ey.Z * (perpx1x3.Z - perpx2x1.Z)
        (D1, D2)

    let compute_jacobians (Dx:Vector<float>) (Dy:Vector<float>) (u:Vector<float>) (v:Vector<float>) =
        let J = CreateMatrix.Dense<float>(2, 2)
        J.[0, 0] <- Dx.DotProduct(u)
        J.[0, 1] <- Dy.DotProduct(u)
        J.[1, 0] <- Dx.DotProduct(v)
        J.[1, 1] <- Dy.DotProduct(v)
        J

    let jacobian_energy A (J:Matrix<float>) =
        let svd_res = svd J
        let sing = svd_res.S
        let s0, s1 = sing.[0], sing.[1]
        let s0sq = s0 * s0
        let s1sq = s1 * s1
        A * (s0sq + 1. / s0sq + s1sq + 1. / s1sq)

    let get_closest_transform (J:Matrix<float>) =
        let svd_res = svd J
        let sing = svd_res.S
        //printfn "J: %A" J
        //printfn "s: %A" sing
        //let s0, s1 = sing.[0], sing.[1]
        //let s0cub = s0 * s0 * s0
        //let s1cub = s1 * s1 * s1
        //let s0_g = 2. * (s0 - 1. / s0cub)
        //let s1_g = 2. * (s1 - 1. / s1cub)
        //let news0 = sqrt (s0_g / (2. * (s0 - 1.)))
        //let news1 = sqrt (s1_g / (2. * (s1 - 1.)))
        let diags = CreateMatrix.DenseOfDiagonalVector(CreateVector.DenseOfArray[|1.; 1.|])
        let MatW = CreateMatrix.DenseIdentity(2)//svd_res.U * diags * (svd_res.U.Transpose ())
        let R = svd_res.U * svd_res.VT
        MatW, R

    let triangle_area = CreateMatrix.DenseOfDiagonalArray [| for i in 0..3 do for t in triangles do yield t.getArea points|]

    let buildA (Dx:Matrix<float>) (Dy:Matrix<float>) (Ws:Matrix<float> array) =
        let IJV = CreateMatrix.Dense<float>(2 * 2 * triangles.Length, 2 * points.Length)
        for row in 0..Dx.RowCount - 1 do
            for col in 0..Dx.ColumnCount - 1 do
                let v = Dx.[row, col]
                if v <> 0. then
                    let W = Ws.[row]
                    IJV.[row, col] <- v * W.[0, 0]
                    IJV.[row, col + points.Length] <- v* W.[0, 1]
                    IJV.[row + 2 * triangles.Length, col] <- v * W.[1, 0]
                    IJV.[row + 2 * triangles.Length, col +  points.Length] <- v * W.[1, 1]
        for row in 0..Dy.RowCount - 1 do
            for col in 0..Dy.ColumnCount - 1 do
                let v = Dx.[row, col]
                let W = Ws.[row]
                if v <> 0. then
                    IJV.[row + triangles.Length, col] <- v * W.[0, 0]
                    IJV.[row + triangles.Length, col + points.Length] <- v * W.[0, 1]
                    IJV.[row + 3 * triangles.Length, col] <- v * W.[1, 0]
                    IJV.[row + 3 * triangles.Length, col +  points.Length] <- v * W.[1, 1]
        IJV

    let build_linear_system (Dx:Matrix<float>) (Dy:Matrix<float>) (W_and_R:(Matrix<float> * Matrix<float>) array)=
        let A = buildA Dx Dy (Array.map fst W_and_R)
        let rhs =
            let trilength = triangles.Length
            let rhs = CreateVector.Dense (2 * 2 * trilength)
            for tri in 0..trilength - 1 do
                let W, R = W_and_R.[tri]
                let tmp = W * R
                rhs.[tri] <- tmp.[0, 0]
                rhs.[tri + trilength] <- tmp.[0, 1]
                rhs.[tri + 2 * trilength] <- tmp.[1, 0]
                rhs.[tri + 3 * trilength] <- tmp.[1, 1]
            rhs
        (A.Transpose () * triangle_area * A), (A.Transpose () * triangle_area * rhs)

    let Dx = CreateMatrix.Dense<float>(triangles.Length, points.Length)
    let Dy = CreateMatrix.Dense<float>(triangles.Length, points.Length)
    triangles |> Array.iteri (fun i t -> let D1, D2 = compute_surface_gradient_matrix t in Dx.SetRow(i, D1); Dy.SetRow(i, D2))

    printfn "Dx:%A" Dx
    printfn "Dy:%A" Dy

    let initial_guess = LSCM points border_point triangles

    let current_u = CreateVector.Dense(initial_guess |> Array.map (function v -> v.X))
    let current_v = CreateVector.Dense(initial_guess |> Array.map (function v -> v.Y))

    let iterations =
        let W_and_R = 
            [| for i in 0..triangles.Length - 1 ->
                let J = compute_jacobians (Dx.Row(i)) (Dy.Row(i)) current_u current_v
                let res = get_closest_transform J
                printfn "W and R %A" res
                res|]
        let A, rhs = build_linear_system Dx Dy W_and_R

        for kv in border_point do
            A.[kv.Key, kv.Key] <- A.[kv.Key, kv.Key] + 1.
            A.[kv.Key + points.Length, kv.Key + points.Length] <- A.[kv.Key + points.Length, kv.Key + points.Length] + 1.
            rhs.[kv.Key] <- rhs.[kv.Key] + kv.Value.X
            rhs.[kv.Key + points.Length] <- rhs.[kv.Key + points.Length] + kv.Value.Y

        printfn "A:%A" A
        printfn "RHS: %A" rhs
        let res = A.Solve(rhs)
        let us = res.[..points.Length - 1]
        let vs = res.[points.Length - 1..]
        let reshaped = [| for i in 0..points.Length - 1 -> Vector2D(us.[i], vs.[i]) |]
        printfn "%A" reshaped
        reshaped

    iterations

