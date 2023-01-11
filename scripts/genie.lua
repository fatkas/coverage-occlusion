MODULE_DIR = path.getabsolute("../")

solution "bgfx"
	configurations {
		"Debug",
		"Release",
	}

	if _ACTION ~= nil and _ACTION:match "^xcode" then
		platforms {
			"Native", -- let xcode decide based on the target output
		}
	else
		platforms {
			"x32",
			"x64",
--			"Xbox360",
			"Native", -- for targets where bitness is not specified
		}
	end

	language "C++"
	startproject "example-occlusion-demo"

BGFX_DIR   = path.getabsolute("../../bgfx")
BX_DIR     = os.getenv("BX_DIR")
BIMG_DIR   = os.getenv("BIMG_DIR")

local BGFX_BUILD_DIR = path.join(MODULE_DIR, "../build")
-- print("build dir " + BGFX_BUILD_DIR)
local BGFX_THIRD_PARTY_DIR = path.join(BGFX_DIR, "3rdparty")
if not BX_DIR then
	BX_DIR = path.getabsolute(path.join(BGFX_DIR, "../bx"))
end

if not BIMG_DIR then
	BIMG_DIR = path.getabsolute(path.join(BGFX_DIR, "../bimg"))
end

if not os.isdir(BX_DIR) or not os.isdir(BIMG_DIR) then

	if not os.isdir(BX_DIR) then
		print("bx not found at \"" .. BX_DIR .. "\". git clone https://github.com/bkaradzic/bx?")
	end

	if not os.isdir(BIMG_DIR) then
		print("bimg not found at \"" .. BIMG_DIR .. "\". git clone https://github.com/bkaradzic/bimg?")
	end

	print("For more info see: https://bkaradzic.github.io/bgfx/build.html")
	os.exit()
end

dofile (path.join(BX_DIR, "scripts/toolchain.lua"))
if not toolchain(BGFX_BUILD_DIR, BGFX_THIRD_PARTY_DIR) then
    print("no toolchain")
	return -- no action specified
end

function copyLib()
end

function exampleProjectDefaults()

	debugdir (path.join(BGFX_DIR, "examples/runtime"))

	includedirs {
		path.join(BIMG_DIR, "include"),
		path.join(BGFX_DIR, "include"),
		path.join(BGFX_DIR, "3rdparty"),
		path.join(BGFX_DIR, "examples/common"),
		path.join(BGFX_DIR, "../enkiTS/src/")
	}

	flags {
		"FatalWarnings",
	}

	links {
		"example-glue",
		"example-common",
		"bgfx",
		"bimg_decode",
		"bimg",
	}

	using_bx()

	configuration { "vs*", "x32 or x64" }
		linkoptions {
			"/ignore:4199", -- LNK4199: /DELAYLOAD:*.dll ignored; no imports found from *.dll
		}
		links { -- this is needed only for testing with GLES2/3 on Windows with VS2008
			"DelayImp",
		}

	configuration { "vs201*", "x32 or x64" }
		linkoptions { -- this is needed only for testing with GLES2/3 on Windows with VS201x
			"/DELAYLOAD:\"libEGL.dll\"",
			"/DELAYLOAD:\"libGLESv2.dll\"",
		}

	configuration { "mingw*" }
		targetextension ".exe"
		links {
			"comdlg32",
			"gdi32",
			"psapi",
		}

	configuration { "vs20*", "x32 or x64" }
		links {
			"gdi32",
			"psapi",
		}

	configuration { "durango" }
		links {
			"d3d11_x",
			"d3d12_x",
			"combase",
			"kernelx",
		}

	configuration { "winstore*" }
		removelinks {
			"DelayImp",
			"gdi32",
			"psapi"
		}
		links {
			"d3d11",
			"d3d12",
			"dxgi"
		}
		linkoptions {
			"/ignore:4264" -- LNK4264: archiving object file compiled with /ZW into a static library; note that when authoring Windows Runtime types it is not recommended to link with a static library that contains Windows Runtime metadata
		}

	-- WinRT targets need their own output directories or build files stomp over each other
	configuration { "x32", "winstore*" }
		targetdir (path.join(BGFX_BUILD_DIR, "win32_" .. _ACTION, "bin", _name))
		objdir (path.join(BGFX_BUILD_DIR, "win32_" .. _ACTION, "obj", _name))

	configuration { "x64", "winstore*" }
		targetdir (path.join(BGFX_BUILD_DIR, "win64_" .. _ACTION, "bin", _name))
		objdir (path.join(BGFX_BUILD_DIR, "win64_" .. _ACTION, "obj", _name))

	configuration { "ARM", "winstore*" }
		targetdir (path.join(BGFX_BUILD_DIR, "arm_" .. _ACTION, "bin", _name))
		objdir (path.join(BGFX_BUILD_DIR, "arm_" .. _ACTION, "obj", _name))

	configuration { "mingw-clang" }
		kind "ConsoleApp"

	configuration { "android*" }
		kind "ConsoleApp"
		targetextension ".so"
		linkoptions {
			"-shared",
		}
		links {
			"EGL",
			"GLESv2",
		}

	configuration { "wasm*" }
		kind "ConsoleApp"

		linkoptions {
			"-s TOTAL_MEMORY=32MB",
			"-s ALLOW_MEMORY_GROWTH=1",
			"--preload-file ../../../examples/runtime@/"
		}

		removeflags {
			"OptimizeSpeed",
		}

		flags {
			"Optimize"
		}

	configuration { "linux-* or freebsd" }
		links {
			"X11",
			"GL",
			"pthread",
		}

	configuration { "rpi" }
		links {
			"X11",
			"brcmGLESv2",
			"brcmEGL",
			"bcm_host",
			"vcos",
			"vchiq_arm",
			"pthread",
		}

	configuration { "osx*" }
		targetdir (path.join(BGFX_BUILD_DIR, "osx_" .. _ACTION))
		objdir (path.join(BGFX_BUILD_DIR, "osx_" .. _ACTION))
		linkoptions {
			"-framework Cocoa",
			"-framework IOKit",
			"-framework OpenGL",
			"-framework QuartzCore",
			"-weak_framework Metal",
		}
		buildoptions_cpp {
			"-Wno-unused-variable",
			"-Wno-shorten-64-to-32",
			"-Wno-unused-function",
		}

	configuration { "ios* or tvos*" }
		targetdir (path.join(BGFX_BUILD_DIR, "ios_" .. _ACTION))
		objdir (path.join(BGFX_BUILD_DIR, "ios_" .. _ACTION))
		kind "ConsoleApp"
		linkoptions {
			"-framework CoreFoundation",
			"-framework Foundation",
			"-framework IOKit",
			"-framework OpenGLES",
			"-framework QuartzCore",
			"-framework UIKit",
			"-weak_framework Metal",
		}

	configuration { "xcode*", "ios" }
		kind "WindowedApp"
		files {
			path.join(BGFX_DIR, "examples/runtime/iOS-Info.plist"),
		}

	configuration { "xcode*", "tvos" }
		kind "WindowedApp"
		files {
			path.join(BGFX_DIR, "examples/runtime/tvOS-Info.plist"),
		}

	configuration {}

	strip()
end

function exampleProject(name)

	project ("example-" .. name)
		uuid (os.uuid("example-" .. name))
		kind "WindowedApp"

	files {
		path.join(MODULE_DIR, "**.c"),
		path.join(MODULE_DIR, "**.cpp"),
		path.join(MODULE_DIR, "**.h"),
		path.join(MODULE_DIR, "../enkiTS/src/TaskScheduler.cpp")
	}

	removefiles {
		path.join(MODULE_DIR, "**.bin.h"),
	}

	defines {
		"ENTRY_CONFIG_IMPLEMENT_MAIN=1",
	}

	exampleProjectDefaults()
end

group "libs"
dofile(path.join(BX_DIR,   "scripts/bx.lua"))
dofile(path.join(BIMG_DIR, "scripts/bimg.lua"))
dofile(path.join(BIMG_DIR, "scripts/bimg_decode.lua"))
dofile(path.join(BGFX_DIR, "scripts/bgfx.lua"))

local function userdefines()
	local defines = {}
	local BGFX_CONFIG = os.getenv("BGFX_CONFIG")
	if BGFX_CONFIG then
		for def in BGFX_CONFIG:gmatch "[^%s:]+" do
			table.insert(defines, "BGFX_CONFIG_" .. def)
		end
	end

	return defines
end

BGFX_CONFIG = userdefines()

bgfxProject("", "StaticLib", BGFX_CONFIG)

if _OPTIONS["with-shared-lib"] then
	group "libs"
	bgfxProject("-shared-lib", "SharedLib", BGFX_CONFIG)
end

if _OPTIONS["with-tools"] then
	group "libs"
	dofile(path.join(BIMG_DIR, "scripts/bimg_encode.lua"))
end

group "examples"
dofile(path.join(BGFX_DIR, "scripts/example-common.lua"))

group "examples"
exampleProject("occlusion-demo")
