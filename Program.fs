// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

open Emgu.CV
open Emgu.CV.Structure
open System.Drawing
open MathNet.Spatial.Euclidean
open triangle
open arap



[<EntryPoint>]
let main argv = 
    use file = CvInvoke.Imread(@"C:\Users\vlj\Desktop\height_map_norway-height-map-aster-30m.png")

    let size = 5
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

    let pinned_val = dict[0, Vector2D(0., 0.); (size * size) - 1, Vector2D(float(size - 1) * step_x, float(size - 1) * step_y)]
        //dict[
        //for i in 0..(size - 1) do
        //    yield (vertex_to_idx i 0, Vector2D(0., float(i) * step_y))
        //for i in 0..(size - 1) do
        //    yield (vertex_to_idx i (size - 1), Vector2D(step_x * float(size), float(i) * step_y))
        //for i in 1..(size - 2) do
        //    yield (vertex_to_idx 0 i, Vector2D(float(i) * step_x, 0.))
        //for i in 1..(size - 2) do
        //    yield (vertex_to_idx (size - 1) i, Vector2D(float(i) * step_x, step_y * float(size)))
        //]


    let pos2D = arap points pinned_val (Array.ofSeq triangles) |> Array.map (function v -> Point(int(v.X), int(v.Y)))

    Array.iteri (fun i p -> CvInvoke.Circle(file, p, 5, Bgr(float(i) * 0., float(i) * 80., 0.).MCvScalar)) pos2D
    let draw_triangle arr = CvInvoke.Polylines(file, Array.ofSeq arr, true, Bgr(255.0, 0.0, 0.0).MCvScalar)
    let to_pts = Seq.map (function (t:triangle) -> t.GetIndexes () |> Array.map (function x -> pos2D.[x])) triangles
    Seq.iter draw_triangle to_pts


    CvInvoke.Imshow("", file)
    CvInvoke.WaitKey 0
    0 // return an integer exit code
