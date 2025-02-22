// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
#version 430 core

layout(location = 0) out vec4 FragColor;
uniform float scale;

in float vdepth;

void main()
{
    //FragColor = vec4(vdepth.xxx * 1.0f/16.0f, 1.0f);
    FragColor = vec4(vdepth.xxx, 1.0f);
}
