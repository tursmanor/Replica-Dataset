// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
#include <EGL.h>
#include <PTexLib.h>
#include <string>
#include <pangolin/image/image_convert.h>
#include <pangolin/image/image_io.h>
#include <pangolin/image/managed_image.h>

#include <Eigen/Geometry>
#include "MirrorRenderer.h"
#include <chrono>
#include <random>
#include <iterator>
#include <iostream>
#include <fstream>

#include<pangolin/geometry/geometry.h>
#include<pangolin/geometry/glgeometry.h>

using namespace std::chrono;

// Type definition for shader options
typedef enum RenderMode { uv=0, tex, color, normal, matcap, vertex, depth, binary, num_modes } RenderType;

/************
  Helpers
************/
// Load and return shader for a given a mode 
pangolin::GlSlProgram loadMeshShaders(const std::string shadir, RenderType mode) {

  const std::string mode_names[] = {"SHOW_UV", "SHOW_TEXTURE", "SHOW_COLOR", "SHOW_NORMAL", "SHOW_MATCAP", "SHOW_VERTEX", "SHOW_DEPTH", "SHOW_BINARY"};
  
  std::map<std::string,std::string> prog_defines;
  for(int i=0; i < (int)RenderMode::num_modes; ++i) {
    prog_defines[mode_names[i]] = std::to_string((int)mode == i);
  }

  pangolin::GlSlProgram shader = pangolin::GlSlProgram();
  shader.AddShaderFromFile(pangolin::GlSlVertexShader, shadir + "/modelviewer.vert", {prog_defines}, {shadir});
  shader.AddShaderFromFile(pangolin::GlSlFragmentShader, shadir + "/modelviewer.frag", {prog_defines}, {shadir});
  shader.Link();
   
  return shader;
}

// Scale, rotate, and translate the mesh by a set of parameters
void getMeshTransformationMatrices(pangolin::OpenGlRenderState s_cam, float headRot, float scaleFactor, float rotation, float translation[], 
                                   pangolin::OpenGlMatrix* objCamMV, pangolin::OpenGlMatrix* objCamKMV){
  // Add dynamic mesh to render
  const pangolin::OpenGlMatrix camMVGL = s_cam.GetModelViewMatrix();
  const Eigen::Matrix<float,4,4> camMV = pangolin::ToEigen<float>( camMVGL );

  const pangolin::OpenGlMatrix camKGL = s_cam.GetProjectionMatrix();
  const Eigen::Matrix<float,4,4> camK = pangolin::ToEigen<float>( camKGL );

  // Create rotation matrix R. X to flip the head up, Z to shake the head left to right.
  Eigen::Matrix<float,4,4> Rx = Eigen::Matrix<float,4,4>();
  Rx(0,0) = 1;   Rx(0,1) = 0;            Rx(0,2) = 0;              Rx(0,3) = 0;
  Rx(1,0) = 0;   Rx(1,1) = cos(headRot); Rx(1,2) = -sin(headRot);  Rx(1,3) = 0;
  Rx(2,0) = 0;   Rx(2,1) = sin(headRot); Rx(2,2) = cos(headRot);   Rx(2,3) = 0;
  Rx(3,0) = 0;   Rx(3,1) = 0;            Rx(3,2) = 0;              Rx(3,3) = 1;

  Eigen::Matrix<float,4,4> Rz = Eigen::Matrix<float,4,4>();
  Rz(0,0) = cos(rotation);   Rz(0,1) = -sin(rotation);  Rz(0,2) = 0;  Rz(0,3) = 0;
  Rz(1,0) = sin(rotation);   Rz(1,1) = cos(rotation);   Rz(1,2) = 0;  Rz(1,3) = 0;
  Rz(2,0) = 0;               Rz(2,1) = 0;               Rz(2,2) = 1;  Rz(2,3) = 0;
  Rz(3,0) = 0;               Rz(3,1) = 0;               Rz(3,2) = 0;  Rz(3,3) = 1;

  Eigen::Matrix<float,4,4> R = Rz * Rx; // Flipped order from the standard bc of row col ordering biz

  // Using 4x4 [R|T] with the right-handed Z rotation matrix
  Eigen::Matrix<float,4,4> objMV = Eigen::Matrix<float,4,4>();
  objMV(0,0) = R(0,0);   objMV(0,1) = R(0,1);  objMV(0,2) = R(0,2);  objMV(0,3) = translation[0];
  objMV(1,0) = R(1,0);   objMV(1,1) = R(1,1);  objMV(1,2) = R(1,2);  objMV(1,3) = translation[1];
  objMV(2,0) = R(2,0);   objMV(2,1) = R(2,1);  objMV(2,2) = R(2,2);  objMV(2,3) = translation[2];
  objMV(3,0) = R(3,0);   objMV(3,1) = R(3,1);  objMV(3,2) = R(3,2);  objMV(3,3) = 1;
       
  // Add scaling matrix
  Eigen::Matrix<float,4,4> S = Eigen::Matrix<float,4,4>();
  S(0,0) = scaleFactor;   S(0,1) = 0;             S(0,2) = 0;           S(0,3) = 0;
  S(1,0) = 0;             S(1,1) = scaleFactor;   S(1,2) = 0;           S(1,3) = 0;
  S(2,0) = 0;             S(2,1) = 0;             S(2,2) = scaleFactor; S(2,3) = 0;
  S(3,0) = 0;             S(3,1) = 0;             S(3,2) = 0;           S(3,3) = 1;

  Eigen::Matrix<float,4,4> res1 = camMV * S * objMV;
  *objCamMV = res1;
  Eigen::Matrix<float,4,4> res2 = camK * camMV * S * objMV;
  *objCamKMV = res2;

}

// Largely taken from PFM_ReadWrite repo on github
// Assumes grayscale input
void saveDepthAsPFM(pangolin::Image<Eigen::Matrix<uint8_t, 3, 1>> image, const std::string path) {

  std::ofstream imageFile(path.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);

  if(imageFile) {

    int width = image.w;
    int height = image.h;
    int numChannels = 1;    // Only one channel bc depth is grayscale

    // Write header
    char type[3];
    type[0] = 'P';
    type[1] = 'f';
    type[2] = 0x0a;
    imageFile << type[0] << type[1] << type[2];

    // Write width and height
    imageFile << width << " " << height << type[2];

    //Assumes little endian storage and ends with a line return 0x0a; stores the type
    char byteOrder[10];
    byteOrder[0] = '-'; byteOrder[1] = '1'; byteOrder[2] = '.'; byteOrder[3] = '0';
    byteOrder[4] = '0'; byteOrder[5] = '0'; byteOrder[6] = '0'; byteOrder[7] = '0';
    byteOrder[8] = '0'; byteOrder[9] = 0x0a;

    for(int i = 0; i < 10; ++i){
      imageFile << byteOrder[i];
    }

    float* buffer = new float[numChannels];

        for(int i = 0 ; i < height ; ++i) {
            for(int j = 0 ; j < width ; ++j) {
                buffer[0] = image(j,height-1-i)[0];
                imageFile.write((char *) buffer, numChannels*sizeof(float));
            }
        }

        delete[] buffer;

        imageFile.close();

  } else {
    std::cerr << "Could not open file for pfm conversion" << std::endl;
  }
}

/***********
  Main
 **********/
int main(int argc, char* argv[]) {

  auto model_start = high_resolution_clock::now();

  const std::string meshFile(argv[1]);
  const std::string atlasFolder(argv[2]);
  const std::string cameraPathFile(argv[4]);
  ASSERT(pangolin::FileExists(meshFile));
  ASSERT(pangolin::FileExists(atlasFolder));
  ASSERT(pangolin::FileExists(cameraPathFile));

  int width = 1920;
  int height = 1080;

  bool renderDepth = true;
  bool saveParameter = true;
  bool renderBinaryMap = true;

  // Machine specific output paths: room name, output save locations, mesh location
  char roomName[10];  char imageOut[1000];  char depthOut[1000]; char paramOut[1000]; char binaryOut[1000];
  strcpy(roomName, meshFile.substr(7,meshFile.length()-16).c_str());	
  strcpy(imageOut, "/home/eleanor/Replica-Dataset/Output/images-path1/%s-%04zu.jpeg");
  strcpy(depthOut, "/home/eleanor/Replica-Dataset/Output/depths-path1/%s-%04zu.pfm");
  strcpy(paramOut, "/home/eleanor/Replica-Dataset/Output/params-path1.txt");
  strcpy(binaryOut, "/home/eleanor/Replica-Dataset/Output/binary-path1/%s-%04zu.jpeg");
  std::string meshPath =  "/home/eleanor/Replica-Dataset/data/lpshead/head.obj";

  // Setup EGL
  EGLCtx egl;
  egl.PrintInformation();
  glFrontFace(GL_CW);     // Don't draw backfaces

  // Setup a framebuffer for rgb image, depth image, binary image
  pangolin::GlTexture render(width, height);
  pangolin::GlRenderBuffer renderBuffer(width, height);
  pangolin::GlFramebuffer frameBuffer(render, renderBuffer);

  pangolin::GlTexture depthTexture(width, height);
  pangolin::GlFramebuffer depthFrameBuffer(depthTexture, renderBuffer);

  pangolin::GlTexture binaryTexture(width,height);
  pangolin::GlFramebuffer binaryFrameBuffer(binaryTexture, renderBuffer);

  // Set up camera path based on script, where the camera is ordered by pox_x,y,z, then look_x,y,z
  //std::vector<float> initCam = {0.5, 0.5, 0.5, -3, 0.5, -0.4}; 
  std::vector<std::vector<float>> cameraPos;
  std::fstream in(cameraPathFile);
  std::string line;
  int i=0;
  
  while(std::getline(in,line)){
    float value;
    std::vector<float> curVec;
    std::stringstream ss(line);
    while(ss >> value){
      curVec.push_back(value);
    }
    cameraPos.push_back(curVec);
    ++i;
  }
  
  int frames = 0;
  pangolin::OpenGlRenderState s_cam(
    pangolin::ProjectionMatrixRDF_BottomLeft(
      width,
      height,
      width / 2.0f,
      width / 2.0f,
      (width - 1.0f) / 2.0f,
      (height - 1.0f) / 2.0f,
      0.1f,
      100.0f),
      pangolin::ModelViewLookAtRDF(cameraPos[frames][0],cameraPos[frames][1],cameraPos[frames][2],cameraPos[frames][3],cameraPos[frames][4],cameraPos[frames][5], 0, 0, 1));

  // Load mirrors
  std::string surfaceFile;
  std::vector<MirrorSurface> mirrors;
  if (surfaceFile.length()) {
    std::ifstream file(surfaceFile);
    picojson::value json;
    picojson::parse(json, file);
    for (size_t i = 0; i < json.size(); i++) {
      mirrors.emplace_back(json[i]);
    }
    std::cout << "Loaded " << mirrors.size() << " mirrors" << std::endl;
  }

  const std::string shadir = STR(SHADER_DIR);
  MirrorRenderer mirrorRenderer(mirrors, width, height, shadir);

  // Load mesh and textures
  PTexMesh ptexMesh(meshFile, atlasFolder, false);

  // Load face mesh and dynamic movement settings for it
  pangolin::Geometry face = pangolin::LoadGeometry(meshPath);
  pangolin::GlGeometry faceGL = pangolin::ToGlGeometry(face);

  // Define shader options
  RenderType type1 = tex;
  RenderType type2 = depth;
  RenderType type3 = binary;
  pangolin::GlSlProgram shader = loadMeshShaders(shadir, type1);
  pangolin::GlSlProgram shaderDepth = loadMeshShaders(shadir, type2);
  pangolin::GlSlProgram shaderBinary = loadMeshShaders(shadir, type3);

  // Initialize face position in room_1
  // Translate: x,y,z
  // Rotate: angle,x,y,z, with x,y,z in [0,1]
  float translation[] = {-1.7,.1,0.2};
  float rotation = M_PI/4; 
  float headRot = M_PI/2;
  float scaleFactor = 2;
  float deltaT = 0.001;
  float deltaR = 0.005; //in radians

  // Initialize output images
  pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> image(width, height);
  pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> depthImage(width, height);
  pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> binaryImage(width, height);  

  // For rendering videos
  int numSpots = 200;
  for(int j = 0; j < numSpots; j++){

    frameBuffer.Bind();

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, width, height);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_CULL_FACE);

    ptexMesh.SetExposure(0.01);

    // Move face into position
    pangolin::OpenGlMatrix objCamMV;
    pangolin::OpenGlMatrix objCamKMV;
    getMeshTransformationMatrices(s_cam, headRot, scaleFactor, rotation, translation, &objCamMV, &objCamKMV);
    
    // Draw face
    shader.Bind();
    shader.SetUniform("T_cam_norm", objCamMV );
    shader.SetUniform("KT_cw", objCamKMV );
    glDisable(GL_CULL_FACE);  // fix texture being inside out
 	  pangolin::GlDraw(shader, faceGL, NULL);
    shader.Unbind();
    
    // Draw scene!
    ptexMesh.Render(s_cam,Eigen::Vector4f(0.0f, 0.0f, 0.0f, 0.0f));

    glPopAttrib(); //GL_VIEWPORT_BIT
    frameBuffer.Unbind();

    // Download and save
    render.Download(image.ptr, GL_RGB, GL_UNSIGNED_BYTE);

    // Modify to fit your file format/where you saved your data
    char imgFilename[1000];
    snprintf(imgFilename, 1000, imageOut, roomName, j);
 
    std::cout << imgFilename << std::endl;
    pangolin::SaveImage(image.UnsafeReinterpret<uint8_t>(),
                        pangolin::PixelFormatFromString("RGB24"),
                        std::string(imgFilename));
   
    // Depth rendering
    if (renderDepth) {
    
      depthFrameBuffer.Bind();

      glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
      glPushAttrib(GL_VIEWPORT_BIT);
      glViewport(0, 0, width, height);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glDisable(GL_CULL_FACE);

      // Draw face
      shaderDepth.Bind();
      shaderDepth.SetUniform("T_cam_norm", objCamMV );
      shaderDepth.SetUniform("KT_cw", objCamKMV );
      pangolin::GlDraw(shaderDepth, faceGL, NULL);
      shaderDepth.Unbind();

      // Draw scene
      ptexMesh.RenderDepth(s_cam,1.f/0.5f,Eigen::Vector4f(0.0f, 0.0f, 0.0f, 0.0f));

      glPopAttrib(); //GL_VIEWPORT_BIT

      depthFrameBuffer.Unbind();
      
      depthTexture.Download(depthImage.ptr, GL_RGB, GL_UNSIGNED_BYTE);

      char filename[1000];
      snprintf(filename, 1000, depthOut, roomName, j);
      saveDepthAsPFM(depthImage, filename);
      //pangolin::SaveImage(
      //    depthImage.UnsafeReinterpret<uint8_t>(),
      //    pangolin::PixelFormatFromString("RGB24"),
      //    std::string(filename));
    }
 
    // Parameter output  
    if (saveParameter) {
      Eigen::Matrix4d T_camera_world = s_cam.GetModelViewMatrix();
      std::vector<float> camPostoSave;

      // Get extrinsics
      camPostoSave.push_back(cameraPos[frames][0]);
      camPostoSave.push_back(cameraPos[frames][1]);
      camPostoSave.push_back(cameraPos[frames][2]);

      // Get intrinsics
      for(int r = 0; r<4; r++){
        for(int c = 0; c<4; c++){
          camPostoSave.push_back(T_camera_world(r,c));
        }
      }

      // Get principal point 
      camPostoSave.push_back(width/2);
      camPostoSave.push_back(height/2);

      // Open new file and append to a new line for each new frame
      std::ofstream myfile(paramOut,std::ios_base::app);
      
      if(myfile){
        for(size_t e = 0; e < camPostoSave.size(); e++){
          myfile << camPostoSave[e]<<" ";
        }
        myfile << std::endl;
        myfile.close();
      }
     
    }

    // Binary map of the foreground/dynamic meshes-- don't draw the background scene 
    if (renderBinaryMap) {
          
      binaryFrameBuffer.Bind();

      glClearColor(0.f, 0.f, 0.f, 0.f);
      glPushAttrib(GL_VIEWPORT_BIT);
      glViewport(0, 0, width, height);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glDisable(GL_CULL_FACE);

      // Draw face
      shaderBinary.Bind();
      shaderBinary.SetUniform("T_cam_norm", objCamMV );
      shaderBinary.SetUniform("KT_cw", objCamKMV );
      pangolin::GlDraw(shaderBinary, faceGL, NULL);
      shaderBinary.Unbind();

      glPopAttrib(); //GL_VIEWPORT_BIT

      binaryFrameBuffer.Unbind();
      
      binaryTexture.Download(binaryImage.ptr, GL_RGB, GL_UNSIGNED_BYTE);

      char binFilename[1000];
      snprintf(binFilename, 1000, binaryOut, roomName, j);
      pangolin::SaveImage(
          binaryImage.UnsafeReinterpret<uint8_t>(),
          pangolin::PixelFormatFromString("RGB24"),
          std::string(binFilename));
    }

    // Modify translation and rotation of face by a small increment
    translation[0] += deltaT;
    rotation += deltaR;
    frames++;

    // Move camera
    if(j<numSpots-1){
      s_cam.SetModelViewMatrix(pangolin::ModelViewLookAtRDF(cameraPos[frames][0],cameraPos[frames][1],cameraPos[frames][2],
                                                            cameraPos[frames][3],cameraPos[frames][4],cameraPos[frames][5],
                                                            0, 0, 1));
    }

    // Pick a random axis and swap rotation on/off every k frames
    if ((frames % 1500) == 0){
      deltaT *= -1;
    }
    if ((frames % 100) == 0){
      deltaR *= -1;
    }
  }

auto model_stop = high_resolution_clock::now();
auto model_duration = duration_cast<microseconds>(model_stop - model_start);
std::cout << "Time taken rendering the model: " << model_duration.count() << " microseconds" << std::endl;
return 0;
}





