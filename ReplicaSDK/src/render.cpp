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
using namespace std::chrono;

int main(int argc, char* argv[]) {

  auto model_start = high_resolution_clock::now();

  ASSERT(argc == 5 || argc == 6, "Usage: ./ReplicaViewer mesh.ply textures [glass.sur] cameraPositions.txt [y if spherical]");

  const std::string meshFile(argv[1]);
  const std::string atlasFolder(argv[2]);
  ASSERT(pangolin::FileExists(meshFile));
  ASSERT(pangolin::FileExists(atlasFolder));

  bool navCam = true;
  std::string navPositions;
  std::string surfaceFile;

  bool spherical = false;

  if (argc == 5) {
    navPositions = std::string(argv[4]);

    if(navPositions.length()==1){

      spherical = true;
      navPositions = std::string(argv[3]);
      ASSERT(pangolin::FileExists(navPositions));

    }else{

      surfaceFile = std::string(argv[3]);
      ASSERT(pangolin::FileExists(surfaceFile));
      ASSERT(pangolin::FileExists(navPositions));

    }
  }

  if (argc == 6) {
    spherical = true;
    surfaceFile = std::string(argv[3]);
    ASSERT(pangolin::FileExists(surfaceFile));

    navPositions = std::string(argv[4]);
    ASSERT(pangolin::FileExists(navPositions));
  }

  // //load nav positions into array
  // std::vector<std::vector<float>> cameraPos;
  // if(navCam){
  //   std::fstream in(navPositions);
  //   std::string line;
  //   int i=0;
  //   while(std::getline(in,line)){
  //     float value;
  //     std::stringstream ss(line);
  //     cameraPos.push_back(std::vector<float>());
  //
  //     int j = 0;
  //     while(ss>>value){
  //       if(j>0){
  //         if(j==1){
  //           cameraPos[i].push_back(value+0.5);
  //         }
  //         else if(j==2){
  //           cameraPos[i].push_back(value-1.0);
  //         }
  //         else if(j==3){
  //           cameraPos[i].push_back(value-1.2);
  //         }else{
  //           cameraPos[i].push_back(value);
  //         }
  //       }
  //       ++j;
  //     }
  //     ++i;
  //   }
  // }

  // load dataset generating txt file
  std::vector<std::vector<float>> cameraPos;
  if(navCam){
    std::fstream in(navPositions);
    std::string line;
    int i=0;
    while(std::getline(in,line)){
      float value;
      std::stringstream ss(line);
      cameraPos.push_back(std::vector<float>());

      while(ss>>value){
        cameraPos[i].push_back(value);
      }
      ++i;
    }
  }


  int width = 1920;
  int height = 1080;
  if(spherical){
    width = 4096;
    height = 2048;
  }

  bool renderDepth = true;
  bool saveParameter = false;
  float depthScale = 65535.0f * 0.1f;

  // Setup EGL
  EGLCtx egl;

  egl.PrintInformation();

  //Don't draw backfaces
  GLenum frontFace = GL_CW;

  if(spherical){
    glFrontFace(frontFace);
  }
  else{
    frontFace = GL_CCW;
    glFrontFace(frontFace);

  }

  // Setup a framebuffer
  pangolin::GlTexture render(width, height);
  pangolin::GlRenderBuffer renderBuffer(width, height);
  pangolin::GlFramebuffer frameBuffer(render, renderBuffer);

  pangolin::GlTexture depthTexture(width, height);
  // pangolin::GlTexture depthTexture(width, height, GL_R32F, false, 0, GL_RED, GL_FLOAT, 0);
  pangolin::GlFramebuffer depthFrameBuffer(depthTexture, renderBuffer);

  // Setup a camera
  std::vector<float> initCam = {0,0.5,-0.6230950951576233};
  if(navCam){
    initCam = cameraPos[0];

    std::cout<<"Initial camera parameters:"<<initCam[0]<<" "<<initCam[1]<<" "<<initCam[2];
  }
  int cx = rand()%4;
  int cy = rand()%4;
  // int cx = 1;
  // int cy = 1;

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

      //for rendering video of camera path
      // pangolin::OpenGlRenderState s_cam(
      //     pangolin::ProjectionMatrixRDF_BottomLeft(
      //         width,
      //         height,
      //         width / 2.0f,
      //         width / 2.0f,
      //         (width - 1.0f) / 2.0f,
      //         (height - 1.0f) / 2.0f,
      //         0.1f,
      //         100.0f),
      //     pangolin::ModelViewLookAtRDF(initCam[0],initCam[1],initCam[2], initCam[0]+1,initCam[1],initCam[2], 0, 0, 1));
      // pangolin::OpenGlRenderState s_cam(
      //     pangolin::ProjectionMatrixRDF_BottomLeft(
      //         width,
      //         height,
      //         width / 2.0f,
      //         width / 2.0f,
      //         (width - 1.0f) / 2.0f,
      //         (height - 1.0f) / 2.0f,
      //         0.1f,
      //         100.0f),
      //     pangolin::ModelViewLookAtRDF(initCam[0],initCam[1],initCam[2], initCam[3],initCam[4],initCam[5], 0, 0, 1));

      // Start at some origin
      Eigen::Matrix4d T_camera_world = s_cam.GetModelViewMatrix();

      // And rotate by 90 degree for each face of the cubemap
      Eigen::Transform<double,3,Eigen::Affine> t(Eigen::AngleAxis<double>(0.5*M_PI,Eigen::Vector3d::UnitY()));
      Eigen::Transform<double,3,Eigen::Affine> u(Eigen::AngleAxis<double>(0.5*M_PI,Eigen::Vector3d::UnitX()));
      Eigen::Transform<double,3,Eigen::Affine> d(Eigen::AngleAxis<double>(M_PI,Eigen::Vector3d::UnitX()));

      Eigen::Matrix4d R_side=Eigen::Matrix4d::Identity();
      Eigen::Matrix4d R_up=Eigen::Matrix4d::Identity();
      Eigen::Matrix4d R_down=Eigen::Matrix4d::Identity();

      R_side=t.matrix();
      R_up=u.matrix();
      R_down=d.matrix();

      //move to the left for new spot
      Eigen::Matrix4d T_new_old = Eigen::Matrix4d::Identity();
      T_new_old.topRightCorner(3, 1) = Eigen::Vector3d(1, 0, 0);

      Eigen::Matrix4d T_stereo = Eigen::Matrix4d::Identity();
      T_stereo.topRightCorner(3, 1) = Eigen::Vector3d(0.2, 0, 0);

      // // get initial rotation matrix
      // Eigen::Matrix3d rot3_before =  Eigen::Quaterniond(initCam[6],initCam[3],initCam[4],initCam[5]).toRotationMatrix();
      // Eigen::Matrix3d rot3_after =  Eigen::Quaterniond(cameraPos[1][6],cameraPos[1][3],cameraPos[1][4],cameraPos[1][5]).toRotationMatrix();
      // Eigen::Matrix4d rot4_before = Eigen::Matrix4d::Identity();
      // Eigen::Matrix4d rot4_after = Eigen::Matrix4d::Identity();
      // rot4_before.block(0,0,3,3) = rot3_before;
      // rot4_after.block(0,0,3,3) = rot3_after;
      //
      // Eigen::Matrix4d rot4 = rot4_after.inverse() * rot4_after;
      // T_camera_world = rot4 * T_camera_world;
      // s_cam.GetModelViewMatrix() = T_camera_world;

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

      // load mesh and textures
      PTexMesh ptexMesh(meshFile, atlasFolder, spherical);

      pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> image(width, height);
      pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> depthImage(width, height);

      // pangolin::ManagedImage<float> depthImage(width, height);
      pangolin::ManagedImage<uint16_t> depthImageInt(width, height);

      // Render 6 frames for cubemap for each spot
      size_t numSpots = 20;
      if(navCam){
        numSpots = cameraPos.size();
      }

      srand(time(0));

      //For rendering equirect videos
      // numSpots = 700;
      std::vector<std::vector<float>> camPostoSave;
      for(size_t step = 0; step<numSpots ; step++){

        Eigen::Matrix4d old_Cam_to_World = s_cam.GetModelViewMatrix();
        size_t eye = 2;
        float basel = 0.f;
        for(int i =0; i<3; i++){
          switch(i){
            case 0:
              eye = 2;
              break;
            case 1:
              eye = 0;
              basel = 0.05;
              break;
            case 2:
              eye = 0;
              basel = 0.1;
              break;
          }
          std::cout<<"Rendering with baseline"<<basel<<std::endl;
          frameBuffer.Bind();
          glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

          glPushAttrib(GL_VIEWPORT_BIT);
          glViewport(0, 0, width, height);
          glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

          glEnable(GL_CULL_FACE);

          ptexMesh.SetExposure(0.01);
          ptexMesh.SetBaseline(basel);
          ptexMesh.Render(s_cam,Eigen::Vector4f(0.0f, 0.0f, 0.0f, 0.0f),eye);

          glDisable(GL_CULL_FACE);

          glPopAttrib(); //GL_VIEWPORT_BIT
          frameBuffer.Unbind();

          for (size_t i = 0; i < mirrors.size(); i++) {
            MirrorSurface& mirror = mirrors[i];
            // capture reflections
            mirrorRenderer.CaptureReflection(mirror, ptexMesh, s_cam, frontFace);

            frameBuffer.Bind();
            glPushAttrib(GL_VIEWPORT_BIT);
            glViewport(0, 0, width, height);

            // render mirror
            mirrorRenderer.Render(mirror, mirrorRenderer.GetMaskTexture(i), s_cam);

            glPopAttrib(); //GL_VIEWPORT_BIT
            frameBuffer.Unbind();
          }

          // Download and save
          render.Download(image.ptr, GL_RGB, GL_UNSIGNED_BYTE);

          char equirectFilename[1000];
          snprintf(equirectFilename, 1000, "/home/selenaling/Desktop/Replica-Dataset/build/ReplicaSDK/equirectData/test-temp1/%s_%04zu_pos%01d.jpeg",meshFile.substr(0,meshFile.length()-9).c_str(),step,i);

          pangolin::SaveImage(image.UnsafeReinterpret<uint8_t>(),
                              pangolin::PixelFormatFromString("RGB24"),
                              std::string(equirectFilename));

          if (renderDepth && i==2) {
            // render depth
            // depthFrameBuffer.Bind();
            // glPushAttrib(GL_VIEWPORT_BIT);
            // glViewport(0, 0, width, height);
            // glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
            //
            // glEnable(GL_CULL_FACE);
            //
            // ptexMesh.RenderDepth(s_cam, depthScale);
            //
            // glDisable(GL_CULL_FACE);
            //
            // glPopAttrib(); //GL_VIEWPORT_BIT
            // depthFrameBuffer.Unbind();
            //
            // depthTexture.Download(depthImage.ptr, GL_RED, GL_FLOAT);
            //
            // // // convert to 16-bit int
            // // for(size_t i = 0; i < depthImage.Area(); i++)
            // //     depthImageInt[i] = static_cast<uint16_t>(depthImage[i] + 0.5f);
            //
            // char filename[1000];
            // snprintf(filename, 1000, "/home/selenaling/Desktop/Replica-Dataset/build/ReplicaSDK/equirectData/test-temp/depth_%s_%04zu.exr",navPositions.substr(0,navPositions.length()-4).c_str(),frame);
            // // pangolin::SaveImage(
            // //     depthImageInt.UnsafeReinterpret<uint8_t>(),
            // //     pangolin::PixelFormatFromString("GRAY16LE"),
            // //     std::string(filename), true, 34.0f);
            // pangolin::SaveImage(
            //   depthImage.UnsafeReinterpret<uint8_t>(),
            //   pangolin::PixelFormatFromString("GRAY32F"),
            //   std::string(filename), pangolin::ImageFileTypeExr, true, 34.0f);

              // //revised
              depthFrameBuffer.Bind();
              glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

              glPushAttrib(GL_VIEWPORT_BIT);
              glViewport(0, 0, width, height);
              glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

              glEnable(GL_CULL_FACE);

              ptexMesh.RenderDepth(s_cam,1.f/0.5f,Eigen::Vector4f(0.0f, 0.0f, 0.0f, 0.0f),eye);

              glDisable(GL_CULL_FACE);

              glPopAttrib(); //GL_VIEWPORT_BIT
              depthFrameBuffer.Unbind();

              depthTexture.Download(depthImage.ptr, GL_RGB, GL_UNSIGNED_BYTE);

              char filename[1000];
              snprintf(filename, 1000, "/home/selenaling/Desktop/Replica-Dataset/build/ReplicaSDK/equirectData/test-temp1/%s_%04zu_pos%01d.jpeg",meshFile.substr(0,meshFile.length()-9).c_str(),step,3);
              pangolin::SaveImage(
                  depthImage.UnsafeReinterpret<uint8_t>(),
                  pangolin::PixelFormatFromString("RGB24"),
                  std::string(filename));

            }

          if (saveParameter) {
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
            camPostoSave[step].push_back(cameraPos[step][0]);
            camPostoSave[step].push_back(cameraPos[step][1]);
            camPostoSave[step].push_back(cameraPos[step][2]);
            // //save the quaternion
            // camPostoSave[step].push_back(q.x());
            // camPostoSave[step].push_back(q.y());
            // camPostoSave[step].push_back(q.z());
            // camPostoSave[step].push_back(q.w());
            for(int r = 0; r<4; r++){
              for(int c = 0; c<4; c++){
                camPostoSave[step].push_back(T_camera_world(r,c));
              }
            }
            camPostoSave[step].push_back(width/2);
            camPostoSave[step].push_back(width/2);

            char parameter_filename[1000];
            snprintf(parameter_filename, 1000, "/home/selenaling/Desktop/Replica-Dataset/dataset/%s_parameters.txt",meshFile.substr(0,meshFile.length()-4).c_str());

            std::ofstream myfile(parameter_filename);
            if(myfile.is_open()){
              for(int step=0;step<numSpots;step++){
                std::cout<<"hello";
                myfile<<step<<" ";

                for(int e = 0; e < camPostoSave[step].size(); e++){
                  myfile<<camPostoSave[step][e]<<" ";
                }
                // myfile<<camPostoSave[step][0]<<" "<<camPostoSave[step][1]<<" "<<camPostoSave[step][2]<<" "
                // <<camPostoSave[step][3]<<" "<<camPostoSave[step][4]<<" "<<camPostoSave[step][5]<<" "<<camPostoSave[step][6]<<" "
                // <<camPostoSave[step][7]<<" "<<camPostoSave[step][8]<<std::endl;
                myfile<<std::endl;
              }
              myfile.close();
            }


          }

        }

        if(step<numSpots-1){
            //set camera to next location in txt file
            s_cam.SetModelViewMatrix(pangolin::ModelViewLookAtRDF(cameraPos[step+1][0],cameraPos[step+1][1],cameraPos[step+1][2],
                                                                  cameraPos[step+1][3],cameraPos[step+1][4],cameraPos[step+1][5],
                                                                  0, 0, 1));
        }


        // For rendering perspective videos
        // size_t numSteps = 320;
        // for(size_t step = 0; step<numSteps ; step++){
        //
        //     Eigen::Matrix4d old_Cam_to_World = s_cam.GetModelViewMatrix();
        //     for(size_t pos = 0; pos<2; pos++){
        //
        //       frameBuffer.Bind();
        //       glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        //
        //       glPushAttrib(GL_VIEWPORT_BIT);
        //       glViewport(0, 0, width, height);
        //       glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        //
        //       glEnable(GL_CULL_FACE);
        //
        //       ptexMesh.SetExposure(0.01);
        //       ptexMesh.Render(s_cam,Eigen::Vector4f(0.0f, 0.0f, 0.0f, 0.0f));
        //
        //       glDisable(GL_CULL_FACE);
        //
        //       glPopAttrib(); //GL_VIEWPORT_BIT
        //       frameBuffer.Unbind();
        //
        //       for (size_t i = 0; i < mirrors.size(); i++) {
        //         MirrorSurface& mirror = mirrors[i];
        //         // capture reflections
        //         mirrorRenderer.CaptureReflection(mirror, ptexMesh, s_cam, frontFace);
        //
        //         frameBuffer.Bind();
        //         glPushAttrib(GL_VIEWPORT_BIT);
        //         glViewport(0, 0, width, height);
        //
        //         // render mirror
        //         mirrorRenderer.Render(mirror, mirrorRenderer.GetMaskTexture(i), s_cam);
        //
        //         glPopAttrib(); //GL_VIEWPORT_BIT
        //         frameBuffer.Unbind();
        //       }
        //
        //       // Download and save
        //       render.Download(image.ptr, GL_RGB, GL_UNSIGNED_BYTE);
        //
        //       char equirectFilename[1000];
        //       snprintf(equirectFilename, 1000, "/home/selenaling/Desktop/Replica-Dataset/build/ReplicaSDK/cubemapData/right_%s_%03zu_pos%01zu.png",navPositions.substr(0,navPositions.length()-4).c_str(),step,pos);
        //
        //       pangolin::SaveImage(
        //           image.UnsafeReinterpret<uint8_t>(),
        //           pangolin::PixelFormatFromString("RGB24"),
        //           std::string(equirectFilename));
        //
        //       if(pos==0){
        //         Eigen::Matrix4d T_target = Eigen::Matrix4d::Identity();
        //         T_target.topRightCorner(3, 1) = Eigen::Vector3d(0.064, 0, 0);
        //
        //         T_camera_world =  T_target.inverse() * old_Cam_to_World;
        //         s_cam.GetModelViewMatrix() = T_camera_world;
        //       }
        //
        //     }
        //
        //     // moving in a square
        //     // s_cam.GetModelViewMatrix() = old_Cam_to_World;
        //     // Eigen::Matrix4d T_target = Eigen::Matrix4d::Identity();
        //     // if(step<numSteps/4){
        //     //   T_target.topRightCorner(3, 1) = Eigen::Vector3d(-0.01, 0, 0);
        //     // }
        //     // else if(step<numSteps/2){
        //     //   T_target.topRightCorner(3, 1) = Eigen::Vector3d(0, 0, 0.01);
        //     // }
        //     // else if(step<numSteps*3/4){
        //     //   T_target.topRightCorner(3, 1) = Eigen::Vector3d(0.01, 0, 0);
        //     // }
        //     // else if(step<numSteps){
        //     //   T_target.topRightCorner(3, 1) = Eigen::Vector3d(0, 0, -0.01);
        //     // }
        //     //
        //     // T_camera_world = T_target.inverse() * old_Cam_to_World ;
        //     // s_cam.GetModelViewMatrix() = T_camera_world;
        //
        //     // moving by the txt file
        //     s_cam.SetModelViewMatrix(pangolin::ModelViewLookAtRDF(cameraPos[step+1][0],cameraPos[step+1][1],cameraPos[step+1][2]+1, 1, 1, cameraPos[step+1][2]+1, 0, 0, 1));
        //
        //
        // }


      std::cout << "\rRendering spot " << numSpots << "/" << numSpots << "... done" << std::endl;
      auto model_stop = high_resolution_clock::now();
      auto model_duration = duration_cast<microseconds>(model_stop - model_start);

      std::cout << "Time taken rendering the model: " << model_duration.count() << " microseconds" << std::endl;
    }
      return 0;
}
