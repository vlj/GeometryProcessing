module arap

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Factorization
open MathNet.Numerics.LinearAlgebra.Matrix
open triangle
open MathNet.Spatial.Euclidean
open System.Numerics
open System.Collections.Generic
open MathNet.Numerics

let closest_rotation A =
    let v = svd A
    let u = v.U
    let vt = v.VT
    if (u * vt).Determinant () < 0. then
        u.SetColumn(u.ColumnCount - 1, -1. * u.Column(u.ColumnCount - 1))
        (u * vt, - v.VT.Transpose () * v.W * v.VT)
    else
        (u * vt, v.VT.Transpose () * v.W * v.VT)

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
    let initial_guess = LSCM points border_point triangles

    let half_edges = dict[
        for t in triangles do
            let [|i0; i1; i2|] = t.GetIndexes ()
            let cotan = t.GetCotangent points
            let [|x0; x1; x2|] = t.project points
            let x0x1 = (x1 - x0).ToVector ()
            let x1x2 = (x2 - x1).ToVector ()
            let x2x0 = (x0 - x2).ToVector ()
            yield ((i0, i1), (x0x1, cotan.[2] / 2.)) // cotan at i2
            yield ((i1, i2), (x1x2, cotan.[0] / 2.)) // cotan at i0
            yield ((i2, i0), (x2x0, cotan.[1] / 2.)) // cotan at i1
        ]

    let insert_or_sum (n: 'a array) k v = n.[k] <- v :: n.[k]

    let source =
        let n = [| for p in points -> ([]: (int * Vector<float> * float) list)|]
        for kv in half_edges do
            let (s, t) = kv.Key
            let (he, cot) = kv.Value
            insert_or_sum n t (s, he, cot)
        n

    let target =
        let n = [| for p in points -> ([]: (int * Vector<float> * float) list)|]
        for kv in half_edges do
            let (s, t) = kv.Key
            let (he, cot) = kv.Value
            insert_or_sum n s (t, he, cot)
        n

    target |>
        Array.iteri (fun i lst->
            lst |> List.iter (function (j, _, _) -> let assertion = List.exists (function (id, _, _) -> id = i) source.[j] in assert assertion)
            )

    let map, pinned_vector = mapping points border_point
    let start_of_fixed = points.Length - border_point.Count

    let pinned_vector2 = CreateMatrix.Dense(pinned_vector.Count, 2)
    for i in 0..pinned_vector.Count - 1 do
        pinned_vector2.[i, 0] <- pinned_vector.[i].Real
        pinned_vector2.[i, 1] <- pinned_vector.[i].Imaginary

    let MatrixA =
        let L = CreateMatrix.Dense<float> (points.Length, points.Length)
        for i in 0..L.RowCount - 1 do
            let assign (j, _, _) =
                let c_l_ij = snd half_edges.[(i, j)]
                if not(half_edges.ContainsKey((j, i))) then
                    L.[i, j] <- - c_l_ij
                else
                    L.[i, j] <- - c_l_ij - snd half_edges.[(j, i)]
            if map.[i] >= start_of_fixed then
                L.[i, i] <- 1.
            else
                target.[i] |> List.iter assign
                let mutable acc = 0.
                for (idx, w) in L.Row(i).EnumerateIndexed () do
                    if idx <> i then
                        acc <- acc + w
                L.[i, i] <- - acc
        L

    let edge_to_triangle =
        dict[
            for t in triangles do
                let [|i0; i1; i2|] = t.GetIndexes ()
                yield ((i0, i1), t)
                yield ((i1, i2), t)
                yield ((i2, i0), t)
            ]

    let find_optimal_Lt (current_uvs: Vector2D array) =
        dict[
            let half_edge_Cs (s, t) =
                let uiuj = (current_uvs.[s] - current_uvs.[t]).ToVector ()
                let (xixj, w) = half_edges.[(s, t)]
                let current_C1 = w * xixj.DotProduct (xixj)
                let current_C2 = w * uiuj.DotProduct(xixj)
                let current_C3 = w * (uiuj.[0] * xixj.[1] - uiuj.[1] * xixj.[0])
                (current_C1, current_C2, current_C3)

            for t in triangles do
                let [|i0; i1; i2|] = t.GetIndexes ()
                let add_half_edge_Cs (C1, C2, C3)(current_C1, current_C2, current_C3) =
                    (C1 + current_C1, C2 + current_C2, C3 + current_C3)
                let (C1, C2, C3) = [(i0, i1); (i1, i2); (i2, i0)] |> List.map half_edge_Cs |> List.fold add_half_edge_Cs (0., 0., 0.)
                let a = C2 / C1
                let b = C3 / C1
                let M = CreateMatrix.DenseOfColumnMajor(2, 2, [a; -b; b; a])
                yield (t, M)
        ]


    let error (current_uvs:Vector2D array) (rotations:IDictionary<triangle, Matrix<float>>) =
        let mutable res = 0.
        let half_edge_energy (s, t) =
            let uiuj = (current_uvs.[s] - current_uvs.[t]).ToVector ()
            let (xixj, w) = half_edges.[(s, t)]
            let rotation = rotations.[edge_to_triangle.[(s, t)]]
            let transformed_xixj = rotation.Multiply(xixj)
            let local_error = (uiuj - transformed_xixj).L2Norm ()
            w * local_error * local_error
        for i in 0..points.Length - 1 do
            let add s (j, (xixj:Vector<float>), w) =
                s + (half_edge_energy (i, j))
            let local_error = (target.[i] |> List.fold add 0.)
            res <- res + local_error
        res

    let iteration current_uvs =
        let rotations = find_optimal_Lt current_uvs

        printfn "before rhs %A" (error current_uvs rotations)

        let rhs =
            let r = CreateMatrix.Dense(points.Length, 2)
            let get_half_bij (r:Matrix<float>) (half_edge:Vector<float>) w = 
                r.Multiply(w / 1.).Multiply(half_edge)
            let get_bij i j =
                let t = let (half_edge, w) = half_edges.[(i, j)] in get_half_bij (rotations.[edge_to_triangle.[(i,j)]]) half_edge w
                if half_edges.ContainsKey (j,i) then
                    let s = let (half_edge, w) = half_edges.[(j, i)] in get_half_bij (rotations.[edge_to_triangle.[(j, i)]]) half_edge w
                    - s + t
                else
                    t
            for i in 0..points.Length - 1 do
                if map.[i] >= start_of_fixed then
                    r.SetRow(i, pinned_vector2.Row(map.[i] - start_of_fixed))
                else
                    let rowcontent = target.[i] |> List.fold (fun s (j, _, _) -> s + (get_bij i j)) (CreateVector.Dense<float> 2)
                    r.SetRow(i, rowcontent)
            r

        let res = MatrixA.Solve(rhs)
        let reshaped = [| for row in res.ToRowArrays () -> Vector2D(row.[0], row.[1]) |]

        printfn "after rhs %A" (error reshaped rotations)
        reshaped

    Seq.fold (fun s _ -> iteration s) initial_guess {0..1}
