// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
#version 430 core

layout(location = 0) out vec4 FragColor;
uniform float scale;

in float depth;

void main()
{
    FragColor = vec4(depth.xxx * 1/16.0f, 1.0f);

}
