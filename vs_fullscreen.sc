$input a_position, a_texcoord0
$output v_texcoord0

/*
 * Copyright 2011-2022 Branimir Karadzic. All rights reserved.
 * License: https://github.com/bkaradzic/bgfx/blob/master/LICENSE
 */

#include "../bgfx/examples/common/common.sh"

void main()
{
	gl_Position = vec4(a_position, 1.0);
	v_texcoord0.x = a_position.x;
	v_texcoord0.y = -a_position.y;
}
