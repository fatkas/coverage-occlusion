## Coverage occlusion culling implementation

it currently only compiles on mac (checkout, run build.sh, it'll fetch dependencies in a upper level folder and also generate xcode solution in build folder)

- working on MT submit/sort/rasterize (sort/rasterize is fine, push is somehow slower)
- working on optimized push function for boxes

### TODO

- make it compilable on windows
- wrap all simd instructions, add neon implementation
