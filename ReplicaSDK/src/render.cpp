// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
#include <EGL.h>
#include <PTexLib.h>
#include <string>
#include <pangolin/image/image_convert.h>
#include <pangolin/image/image_io.h>

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

int main(int argc, char* argv[]) {

  auto model_start = high_resolution_clock::now();

  const std::string meshFile(argv[1]);
  const std::string atlasFolder(argv[2]);
  ASSERT(pangolin::FileExists(meshFile));
  ASSERT(pangolin::FileExists(atlasFolder));

  std::string surfaceFile;

  int width = 1920;
  int height = 1080;

  bool renderDepth = true;
  bool saveParameter = false;
  float depthScale = 65535.0f * 0.1f;

  // Setup EGL
  EGLCtx egl;
  egl.PrintInformation();
  glFrontFace(GL_CW);   //Don't draw backfaces

  // Setup a framebuffer
  pangolin::GlTexture render(width, height);
  pangolin::GlRenderBuffer renderBuffer(width, height);
  pangolin::GlFramebuffer frameBuffer(render, renderBuffer);

  pangolin::GlTexture depthTexture(width, height);
  pangolin::GlFramebuffer depthFrameBuffer(depthTexture, renderBuffer);

  // Set up a camera
  // for room_1
  std::vector<float> initCam = {0.5, 0.5, 0.5, -3, 0.5, -0.4}; // pos_x,y,z then look_x,y,z

  //for room_2
  // std::vector<float> initCam = {0.8409916162490845,0.605334997177124,-1.0, 1, 1, -1.0};

  // //for office_0
  // std::vector<float> initCam = {-1.2106552124023438,-0.4804558753967285,0.38081249594688416, 1, 1, 0.38081249594688416 };

  //for apartment_0
  //std::vector<float> initCam = {-2.4934263229370117,-4.147322177886963,-0.0, 1, 1, -0.0};
  
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
      pangolin::ModelViewLookAtRDF(initCam[0],initCam[1],initCam[2], initCam[3],initCam[4], initCam[5], 0, 0, 1));

  // Start at some origin
  Eigen::Matrix4d T_camera_world = s_cam.GetModelViewMatrix();

  // And rotate by 90 degree for each face of the cubemap for cubemap rendering
  Eigen::Transform<double,3,Eigen::Affine> t(Eigen::AngleAxis<double>(0.5*M_PI,Eigen::Vector3d::UnitY()));
  Eigen::Transform<double,3,Eigen::Affine> u(Eigen::AngleAxis<double>(0.5*M_PI,Eigen::Vector3d::UnitX()));
  Eigen::Transform<double,3,Eigen::Affine> d(Eigen::AngleAxis<double>(M_PI,Eigen::Vector3d::UnitX()));

  Eigen::Matrix4d R_side=Eigen::Matrix4d::Identity();
  Eigen::Matrix4d R_up=Eigen::Matrix4d::Identity();
  Eigen::Matrix4d R_down=Eigen::Matrix4d::Identity();

  R_side = t.matrix();
  R_up = u.matrix();
  R_down = d.matrix();

  // load mirrors
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
  pangolin::Geometry face = pangolin::LoadGeometry("/home/eleanor/Replica-Dataset/data/lpshead/head.obj");
  pangolin::GlGeometry faceGL = pangolin::ToGlGeometry(face);

  // Define shader options
  enum class RenderMode2 { uv=0, tex, color, normal, matcap, vertex, depth, num_modes };
  const std::string mode_names[] = {"SHOW_UV", "SHOW_TEXTURE", "SHOW_COLOR", "SHOW_NORMAL", "SHOW_MATCAP", "SHOW_VERTEX", "SHOW_DEPTH"};
  std::map<std::string,std::string> prog_defines;
  for(int i=0; i < (int)RenderMode2::num_modes; ++i) {
      std::cout << i << std::endl;
      prog_defines[mode_names[i]] = std::to_string((int)RenderMode2::tex == i);
  }

  // Load shader for face
  pangolin::GlSlProgram shader = pangolin::GlSlProgram();
  shader.AddShaderFromFile(pangolin::GlSlVertexShader, shadir + "/modelviewer.vert", {prog_defines}, {shadir});
  shader.AddShaderFromFile(pangolin::GlSlFragmentShader, shadir + "/modelviewer.frag", {prog_defines}, {shadir});
  shader.Link();

  // Define and load depth shader for face
  std::map<std::string,std::string> prog_definesDepth;
  for(int i=0; i < (int)RenderMode2::num_modes; ++i) {
      prog_definesDepth[mode_names[i]] = std::to_string((int)RenderMode2::depth == i);
  }

  // Load shader for face
  pangolin::GlSlProgram shaderDepth = pangolin::GlSlProgram();
  shaderDepth.AddShaderFromFile(pangolin::GlSlVertexShader, shadir + "/modelviewer.vert", {prog_definesDepth}, {shadir});
  shaderDepth.AddShaderFromFile(pangolin::GlSlFragmentShader, shadir + "/modelviewer.frag", {prog_definesDepth}, {shadir});
  shaderDepth.Link();


  // Initialize face position in room_1
  // Translate: x,y,z
  // Rotate: angle,x,y,z, with x,y,z in [0,1]
  float translation[] = {-1.7,.1,0.2};
  float rotation = M_PI/4; 
  float headRot = M_PI/2;
  float scaleFactor = 2;
  float deltaT = 0.001;
  float deltaR = 0.005; //in radians
  int frames = 0;

  // Initialize output images
  pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> image(width, height);
  pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> depthImage(width, height);
  pangolin::ManagedImage<uint16_t> depthImageInt(width, height);

  // Output room name
  //std::string roomName = meshFile.substr(7,meshFile.length()-16).c_str();

  // For rendering equirect videos
  int numSpots = 200;
  std::vector<std::vector<float>> camPostoSave;
  for(size_t j = 0; j<numSpots; j++){

    frameBuffer.Bind();

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glPushAttrib(GL_VIEWPORT_BIT);
    glViewport(0, 0, width, height);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glEnable(GL_CULL_FACE);

    ptexMesh.SetExposure(0.01);
    //ptexMesh.SetBaseline(basel);

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

    const Eigen::Matrix<float,4,4> res1 = camMV * S * objMV;
    const pangolin::OpenGlMatrix objCamMV = pangolin::OpenGlMatrix( res1 );
    const Eigen::Matrix<float,4,4> res2 = camK * camMV * S * objMV;
    const pangolin::OpenGlMatrix objCamKMV = pangolin::OpenGlMatrix( res2 );
 
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
    char equirectFilename[1000];
    snprintf(equirectFilename, 1000, "/home/eleanor/Replica-Dataset/test_video/%s-%04zu.jpeg", meshFile.substr(7,meshFile.length()-16).c_str(), j);
 
    std::cout << equirectFilename << std::endl;
    pangolin::SaveImage(image.UnsafeReinterpret<uint8_t>(),
                        pangolin::PixelFormatFromString("RGB24"),
                        std::string(equirectFilename));
   
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
      snprintf(filename, 1000,"/home/eleanor/Replica-Dataset/test_depth/%s-%04zu.jpeg", meshFile.substr(7,meshFile.length()-16).c_str(), j);
      pangolin::SaveImage(
          depthImage.UnsafeReinterpret<uint8_t>(),
          pangolin::PixelFormatFromString("RGB24"),
          std::string(filename));
      
    }
 
      // Parameter output  
      /* if (saveParameter) {
        //calculate and save the quaternion
        T_camera_world = s_cam.GetModelViewMatrix();
        Eigen::Matrix4d cam_matrix = T_camera_world.inverse();
        Eigen::Vector3d forward_after = cam_matrix.block(0,2,3,1);
        Eigen::Vector3d right_after = cam_matrix.block(0,0,3,1);
        Eigen::Vector3d up_after = cam_matrix.block(0,1,3,1);

        right_after.normalized();
        double theta = forward_after.dot(up_after);
        double angle_rotation = -1.0*acos(theta);

        std::cout<<"quaternion matrix"<<std::endl;
        Eigen::Quaterniond q;
        q = Eigen::AngleAxis<double>(angle_rotation,right_after);
        std::cout<<q.w()<<" "<<q.x()<<" "<<q.y()<<" "<<q.z()<<std::endl;
        std::cout<<T_camera_world<<std::endl;
        for(int r = 0; r<4; r++){
          for(int c = 0; c<4; c++){
            std::cout<<T_camera_world(r,c);
          }
          std::cout<<std::endl;
        }

        camPostoSave.push_back(std::vector<float>());
        camPostoSave[j].push_back(cameraPos[j][0]);
        camPostoSave[j].push_back(cameraPos[j][1]);
        camPostoSave[j].push_back(cameraPos[j][2]);
        // //save the quaternion
        // camPostoSave[step].push_back(q.x());
        // camPostoSave[step].push_back(q.y());
        // camPostoSave[step].push_back(q.z());
        // camPostoSave[step].push_back(q.w());
        for(int r = 0; r<4; r++){
          for(int c = 0; c<4; c++){
            camPostoSave[j].push_back(T_camera_world(r,c));
          }
        }
        camPostoSave[j].push_back(width/2);
        camPostoSave[j].push_back(width/2);

        char parameter_filename[1000];
        snprintf(parameter_filename, 1000, "/home/eleanor/Replica-Dataset/dataset/%s_parameters.txt",meshFile.substr(0,meshFile.length()-4).c_str());

        std::ofstream myfile(parameter_filename);
        if(myfile.is_open()){
          for(int j=0;j<numSpots;j++){
            std::cout<<"hello";
            myfile<<j<<" ";

            for(int e = 0; e < camPostoSave[j].size(); e++){
              myfile<<camPostoSave[j][e]<<" ";
            }
            // myfile<<camPostoSave[step][0]<<" "<<camPostoSave[step][1]<<" "<<camPostoSave[step][2]<<" "
            // <<camPostoSave[step][3]<<" "<<camPostoSave[step][4]<<" "<<camPostoSave[step][5]<<" "<<camPostoSave[step][6]<<" "
            // <<camPostoSave[step][7]<<" "<<camPostoSave[step][8]<<std::endl;
            myfile<<std::endl;
          }
          myfile.close();
        } 

      }*/

  if(j<numSpots-1){
      // Move camera a small amount in x + look up in z
      float deltaMove = 1.f/100.f;
      float deltaLook = 1.f/250.f;
      s_cam.SetModelViewMatrix(pangolin::ModelViewLookAtRDF(initCam[0] -= deltaMove,initCam[1],initCam[2],
                                                            initCam[3],initCam[4],initCam[5] += deltaLook,
                                                            0, 0, 1));
        
    }

    // Modify translation and rotation of face by a small increment
    translation[0] += deltaT;
    rotation += deltaR;
    frames++;

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