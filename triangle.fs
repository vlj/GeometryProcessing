module triangle

open MathNet.Spatial.Euclidean

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
            let BA = - AB
            let CB = - BC
            let AC = - CA
            let a = AB.Length
            let b = BC.Length
            let c = CA.Length
            let args = [| AB, AC, c * b; BA, BC, c * a; CA, CB, a * b |]
            let coss = args |> Array.map (function (e1, e2, denum) -> e1.DotProduct(e2) / denum)
            let sins = args |> Array.map (function (e1, e2, denum) -> e1.CrossProduct(e2).Length / denum)
            Array.map2 (/) coss sins
        member this.project (vertexes: Vector3D array) =
            let (Ex, Ey, _) = this.GetLocalBasis vertexes
            let (p0, p1, p2) = this.GetVertex vertexes
            [| Vector2D(p0.DotProduct(Ex), p0.DotProduct(Ey));
               Vector2D(p1.DotProduct(Ex), p1.DotProduct(Ey));
               Vector2D(p2.DotProduct(Ex), p2.DotProduct(Ey))|]
        member this.getGradient (vertexes: Vector3D array) =
            let [|i1;i2;i3|] = this.GetIndexes ()
            let x3x2 = vertexes.[i3] - vertexes.[i2]
            let x1x3 = vertexes.[i1] - vertexes.[i3]
            let x2x1 = vertexes.[i2] - vertexes.[i1]

            let n = x3x2.CrossProduct(x1x3)
            // twice the triangle area
            let Ax2 = sqrt n.Length
            let u = n.Normalize()

            let perpx2x1 = x2x1.Length * u.CrossProduct(x2x1).Normalize() / Ax2
            let perpx1x3 = x1x3.Length * u.CrossProduct(x1x3).Normalize() / Ax2

            perpx1x3, perpx2x1
