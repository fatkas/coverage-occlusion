cd ..
git clone https://github.com/bkaradzic/bx.git
git clone https://github.com/bkaradzic/bimg.git
git clone https://github.com/bkaradzic/bgfx.git

cd bgfx
make -j8 shaderc

cd ../occlusion-demo
../bx/tools/bin/darwin/genie --xcode=osx xcode11
make rebuild