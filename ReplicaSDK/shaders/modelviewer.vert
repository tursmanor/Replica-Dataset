#version 130

#expect SHOW_COLOR
#expect SHOW_NORMAL
#expect SHOW_TEXTURE
#expect SHOW_MATCAP
#expect SHOW_UV
#expect SHOW_DEPTH

    uniform mat4 T_cam_norm;
    uniform mat4 KT_cw;
    attribute vec3 vertex;
    varying float depth; 
    uniform vec4 clipPlane;
    vec4 position = vec4(vertex, 1.0);

#if SHOW_COLOR
    attribute vec4 color;
    varying vec4 vColor;
    void main() {
        vColor = color;
#elif SHOW_NORMAL
    attribute vec3 normal;
    varying vec3 vNormal;
    void main() {
        vNormal = mat3(T_cam_norm) * normal;
#elif SHOW_TEXTURE
    attribute vec2 uv;
    varying vec2 vUV;
    void main() {
        vUV = vec2(uv.s, 1.0-uv.t); // try flipping to get face texture correct
#elif SHOW_MATCAP
    attribute vec3 normal;
    varying vec3 vNormalCam;
    void main() {
        vNormalCam = mat3(T_cam_norm) * normal;
#elif SHOW_UV
    attribute vec2 uv;
    varying vec2 vUV;
    void main() {
        vUV = uv;
#elif SHOW_DEPTH
    vec4 cameraPos = T_cam_norm * position;
    void main() {
        depth = cameraPos.z;
        gl_ClipDistance[0] = dot(position, clipPlane);
#else
    varying vec3 vP;
    void main() {
        vP = vertex;
#endif

    gl_Position = KT_cw * position;
}