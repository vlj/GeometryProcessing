module slim

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Factorization
open MathNet.Numerics.LinearAlgebra.Matrix
open triangle
open MathNet.Spatial.Euclidean
open System.Numerics
open System.Collections.Generic
open MathNet.Numerics

let SLIM (points : Vector3D array) (border_point: IDictionary<int, Vector2D>) (triangles: triangle array) =
    let compute_surface_gradient_matrix (t:triangle) =
        let [|i1;i2;i3|] = t.GetIndexes ()
        let perpx1x3, perpx2x1 = t.getGradient points
        let (Ex, Ey, _) = t.GetLocalBasis points
        let D1 = CreateVector.Dense<float> (1)
        D1.[i2] <- Ex.X * perpx1x3.X + Ex.Y * perpx1x3.Y + Ex.Z * perpx1x3.Z
        D1.[i3] <- Ex.X * perpx2x1.X + Ex.Y * perpx2x1.Y + Ex.Z * perpx2x1.Z
        D1.[i1] <- - Ex.X * (perpx1x3.X - perpx2x1.X) - Ex.Y * (perpx1x3.Y - perpx2x1.Y) - Ex.Z * (perpx1x3.Z - perpx2x1.Z)
        let D2 = CreateVector.Dense<float> (1)
        D2.[i2] <- Ey.X * perpx1x3.X + Ey.Y * perpx1x3.Y + Ey.Z * perpx1x3.Z
        D2.[i3] <- Ey.X * perpx2x1.X + Ey.Y * perpx2x1.Y + Ey.Z * perpx2x1.Z
        D2.[i1] <- - Ey.X * (perpx1x3.X - perpx2x1.X) - Ey.Y * (perpx1x3.Y - perpx2x1.Y) - Ey.Z * (perpx1x3.Z - perpx2x1.Z)
        (D1, D2)

    let compute_jacobians (Dx:Vector<float>) (Dy:Vector<float>) (u:Vector<float>) (v:Vector<float>) =
        let J = CreateMatrix.Dense<float>(2, 2)
        J.[0, 0] <- Dx.DotProduct(u)
        J.[0, 1] <- Dy.DotProduct(v)
        J.[1, 0] <- Dx.DotProduct(u)
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
        let s0, s1 = sing.[0], sing.[1]
        let s0cub = s0 * s0 * s0
        let s1cub = s1 * s1 * s1
        let s0_g = 2. * (s0 - 1. / s0cub)
        let s1_g = 2. * (s1 - 1. / s1cub)
        let news0 = sqrt (s0_g / (2. * (s0 - 1.)))
        let news1 = sqrt (s1_g / (2. * (s1 - 1.)))
        let diags = CreateMatrix.DenseOfDiagonalVector(CreateVector.DenseOfArray[|news0; news1|])
        let MatW = svd_res.U * diags * (svd_res.U.Transpose ())
        let R = svd_res.U * svd_res.VT
        MatW, R

    let buildA (Dx:Matrix<float>) (Dy:Matrix<float>) (W_11:Vector<float>) (W_12:Vector<float>) (W_21:Vector<float>) (W_22:Vector<float>) =
        let IJV = CreateMatrix.Dense<float>(2 * 2 * triangles.Length, 2 * points.Length)
        for row in 0..Dx.RowCount - 1 do
            for col in 0..Dx.ColumnCount - 1 do
                let v = Dx.[row, col]
                IJV.[row, col] <- v * W_11.[row]
                IJV.[row, col + points.Length] <- v * W_12.[row]
                IJV.[row + 2 * triangles.Length, col] <- v * W_21.[row]
                IJV.[row + 2 * triangles.Length, col +  points.Length] <- v * W_22.[row]
        for row in 0..Dy.RowCount - 1 do
            for col in 0..Dy.ColumnCount - 1 do
                let v = Dx.[row, col]
                IJV.[row, col] <- v * W_11.[row]
                IJV.[row, col + points.Length] <- v * W_12.[row]
                IJV.[row + 2 * triangles.Length, col] <- v * W_21.[row]
                IJV.[row + 2 * triangles.Length, col +  points.Length] <- v * W_22.[row]
        IJV

    let build_linear_system (Dx:Matrix<float>) (Dy:Matrix<float>) (W_11:Vector<float>) (W_12:Vector<float>) (W_21:Vector<float>) (W_22:Vector<float>) =
        let A = buildA Dx Dy W_11 W_12 W_21 W_22
        let buildRhs (R:Matrix<float>) =
            let trilength = triangles.Length
            let rhs = CreateVector.Dense (2 * 2 * trilength)
            for tri in 0..trilength - 1 do
                rhs.[tri] <- W_11.[tri] * R.[tri, 0] + W_12.[tri] * R.[tri, 1]
                rhs.[tri + trilength] <- W_11.[tri] * R.[tri, 2] + W_12.[tri] * R.[tri, 3]
                rhs.[tri + 2 * trilength] <- W_21.[tri] * R.[tri, 0] + W_22.[tri] * R.[tri, 1]
                rhs.[tri + 3 * trilength] <- W_21.[tri] * R.[tri, 2] + W_22.[tri] * R.[tri, 3]

        A.Transpose()

    let Dx = CreateMatrix.Dense<float>(triangles.Length, points.Length)
    let Dy = CreateMatrix.Dense<float>(triangles.Length, points.Length)
    triangles |> Array.iteri (fun i t -> let D1, D2 = compute_surface_gradient_matrix t in Dx.SetRow(i, D1); Dy.SetRow(i, D2))

    let current_u = CreateVector.Dense(points.Length)
    let current_v = CreateVector.Dense(points.Length)

    let iterations =
        for i in 0..triangles.Length - 1 do
            let J = compute_jacobians (Dx.Row(i)) (Dy.Row(i)) current_u current_v
            let tmp = get_closest_transform J
            tmp

        ()

    ()

