// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
#include <PTexLib.h>
#include <ctime>
#include <pangolin/display/display.h>
#include <pangolin/display/widgets/widgets.h>

#include <pangolin/pangolin.h>
#include <pangolin/geometry/geometry_obj.h>

#include <pangolin/geometry/glgeometry.h>
#include <math.h>

#include "MirrorRenderer.h"

int main(int argc, char* argv[]) {
  ASSERT(argc==3 || argc == 4 || argc == 5, "Usage: ./ReplicaViewer mesh.ply textures [glass.sur] [y if spherical]");

  const std::string meshFile(argv[1]);
  const std::string atlasFolder(argv[2]);
  ASSERT(pangolin::FileExists(meshFile));
  ASSERT(pangolin::FileExists(atlasFolder));

  std::string surfaceFile;
  bool spherical = false;

  if (argc == 4) {
    surfaceFile = std::string(argv[3]);
    if(surfaceFile.length()==1){
      spherical=true;
    }else{
      ASSERT(pangolin::FileExists(surfaceFile));

    }
  }

  if (argc==5) {

    surfaceFile = std::string(argv[3]);
    ASSERT(pangolin::FileExists(surfaceFile));

    spherical = true;
  }

  const int uiWidth = 180;
  const int width = 1280;
  const int height = 960;

  // Setup OpenGL Display (based on GLUT)
  pangolin::CreateWindowAndBind("ReplicaViewer", uiWidth + width, height);

  if (glewInit() != GLEW_OK) {
    pango_print_error("Unable to initialize GLEW.");
  }

  // Setup default OpenGL parameters
  glEnable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  const GLenum frontFace = GL_CW;
  glFrontFace(frontFace);
  glLineWidth(1.0f);

  // Tell the base view to arrange its children equally
  if (uiWidth != 0) {
    pangolin::CreatePanel("ui").SetBounds(0, 1.0f, 0, pangolin::Attach::Pix(uiWidth));
  }

  pangolin::View& container =
      pangolin::CreateDisplay().SetBounds(0, 1.0f, pangolin::Attach::Pix(uiWidth), 1.0f);

  pangolin::OpenGlRenderState s_cam(
      pangolin::ProjectionMatrixRDF_TopLeft(
          width,
          height,
          width / 2.0f,
          width / 2.0f,
          (width - 1.0f) / 2.0f,
          (height - 1.0f) / 2.0f,
          0.1f,
          100.0f),
      pangolin::ModelViewLookAtRDF(0, 1, 0.1, 1, 1, 0.1, 0, 0, 1));


  pangolin::Handler3D s_handler(s_cam);

  pangolin::View& meshView = pangolin::Display("MeshView")
                                 .SetBounds(0, 1.0f, 0, 1.0f, (double)width / (double)height)
                                 .SetHandler(&s_handler);

  container.AddDisplay(meshView);

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

  pangolin::Var<float> exposure("ui.Exposure", 0.01, 0.0f, 0.1f);
  pangolin::Var<float> gamma("ui.Gamma", ptexMesh.Gamma(), 1.0f, 3.0f);
  pangolin::Var<float> saturation("ui.Saturation", ptexMesh.Saturation(), 0.0f, 2.0f);
  pangolin::Var<float> depthScale("ui.Depth_scale", 0.1f, 0.0f, 1.0f);

  pangolin::Var<bool> wireframe("ui.Wireframe", false, true);
  pangolin::Var<bool> renderImage("ui.Render_view", false, true);
  pangolin::Var<bool> drawBackfaces("ui.Draw_backfaces", false, true);
  pangolin::Var<bool> drawMirrors("ui.Draw_mirrors", true, true);
  pangolin::Var<bool> drawDepth("ui.Draw_depth", false, true);

  ptexMesh.SetExposure(exposure);

  // Try adding the bunny mesh to the scene
  pangolin::Geometry bunny = pangolin::LoadGeometry("/home/eleanor/Replica-Dataset/data/lpshead/head.obj");
  pangolin::GlGeometry bunnyGL = pangolin::ToGlGeometry(bunny);

  // Define shader options
  enum class RenderMode { uv=0, tex, color, normal, matcap, vertex, num_modes };
  const std::string mode_names[] = {"SHOW_UV", "SHOW_TEXTURE", "SHOW_COLOR", "SHOW_NORMAL", "SHOW_MATCAP", "SHOW_VERTEX"};
  std::map<std::string,std::string> prog_defines;
  for(int i=0; i < (int)RenderMode::num_modes-1; ++i) {
      prog_defines[mode_names[i]] = std::to_string((int)RenderMode::tex == i);
  }

  // Load shader for bunny
  pangolin::GlSlProgram shader = pangolin::GlSlProgram();
  shader.AddShaderFromFile(pangolin::GlSlVertexShader, shadir + "/modelviewer.vert", {prog_defines}, {shadir});
  shader.AddShaderFromFile(pangolin::GlSlFragmentShader, shadir + "/modelviewer.frag", {prog_defines}, {shadir});
  shader.Link();

  // Initialize bunny position
  // Translate: x,y,z
  // Rotate: angle,x,y,z, with x,y,z in [0,1]
  float translation[] = {-2.5,0.5,0};
  float rotation = M_PI/4; 
  float headRot = M_PI/2;
  float deltaT = 0.001;
  float deltaR = 0.005; //in radians
  int frames = 0;

  while (!pangolin::ShouldQuit()) {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    if (exposure.GuiChanged()) {
      ptexMesh.SetExposure(exposure);
    }

    if (gamma.GuiChanged()) {
      ptexMesh.SetGamma(gamma);
    }

    if (saturation.GuiChanged()) {
      ptexMesh.SetSaturation(saturation);
    }

    if (meshView.IsShown()) {
      meshView.Activate(s_cam);

      if(renderImage){

        pangolin::GlTexture render(width, height);
        pangolin::GlRenderBuffer renderBuffer(width, height);
        pangolin::GlFramebuffer frameBuffer(render, renderBuffer);

        pangolin::ManagedImage<Eigen::Matrix<uint8_t, 3, 1>> image(width, height);

        frameBuffer.Bind();
        glPushAttrib(GL_VIEWPORT_BIT);
        glViewport(0, 0, width, height);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

        glEnable(GL_CULL_FACE);

        ptexMesh.Render(s_cam);
        glDisable(GL_CULL_FACE);

        glPopAttrib(); //GL_VIEWPORT_BIT
        frameBuffer.Unbind();
        render.Download(image.ptr, GL_RGB, GL_UNSIGNED_BYTE);

        char cmapFilename[1000];
        snprintf(cmapFilename, 1000, "/home/eleanor/Replica-Dataset/output/renderFromViewer.png");

        pangolin::SaveImage(
            image.UnsafeReinterpret<uint8_t>(),
            pangolin::PixelFormatFromString("RGB24"),
            std::string(cmapFilename));
      }

      if (drawBackfaces) {
        glDisable(GL_CULL_FACE);
      } else {
        glEnable(GL_CULL_FACE);
      }

      if (wireframe) {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        ptexMesh.Render(s_cam);
        glDisable(GL_POLYGON_OFFSET_FILL);
        // render wireframe on top
        ptexMesh.RenderWireframe(s_cam);
      } 
      else if (drawDepth) {
        ptexMesh.RenderDepth(s_cam, depthScale);
      } 
      else {

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
        
        const Eigen::Matrix<float,4,4> res1 = camMV * objMV;
        const pangolin::OpenGlMatrix objCamMV = pangolin::OpenGlMatrix( res1 );
        const Eigen::Matrix<float,4,4> res2 = camK * camMV * objMV;
        const pangolin::OpenGlMatrix objCamKMV = pangolin::OpenGlMatrix( res2 );

        shader.Bind();
        shader.SetUniform("T_cam_norm", objCamMV );
        shader.SetUniform("KT_cw", objCamKMV );
        glDisable(GL_CULL_FACE);  // fix texture being inside out
   	    pangolin::GlDraw(shader, bunnyGL, NULL);
        shader.Unbind();

        // Draw scene!
        ptexMesh.Render(s_cam);
      
        // Modify translation and rotation by a small increment
        translation[0] += deltaT;
        rotation += deltaR;
        frames++;

        // Pick a random axis and swap rotation on/off every k frames
        if ((frames % 2000) == 0){
          deltaT *= -1;
        }
        if ((frames % 100) == 0){
          deltaR *= -1;
        }
      }

      glDisable(GL_CULL_FACE);

      if (drawMirrors) {
        for (size_t i = 0; i < mirrors.size(); i++) {
          MirrorSurface& mirror = mirrors[i];
          // capture reflections
          mirrorRenderer.CaptureReflection(mirror, ptexMesh, s_cam, frontFace, drawDepth, depthScale);

          // render mirror
          mirrorRenderer.Render(mirror, mirrorRenderer.GetMaskTexture(i), s_cam, drawDepth);
        }
      }

    }

    pangolin::FinishFrame();
  }

  return 0;
}
