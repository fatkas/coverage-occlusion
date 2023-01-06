#
# Copyright 2011-2022 Branimir Karadzic. All rights reserved.
# License: https://github.com/bkaradzic/bgfx/blob/master/LICENSE
#

TEMP=../build

include ../bgfx/scripts/shader-embeded.mk

rebuild:
	@make -s --no-print-directory clean all
