#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iostream>
#define STB_IMAGE_IMPLEMENTATION
const double PI = 3.141592653;

#include "box.h"
#include "bvh.h"
#include "camera.h"
#include "hitablelist.h"
#include "jyorand.h"
#include "kuinkerm.h"
#include "material.h"
#include "onb.h"
#include "pdf.h"
#include "perlin.h"
#include "rectangle.h"
#include "smoke.h"
#include "sphere.h"
#include "stb_image.h"
#include "texture.h"
#include "transformation.h"
#include "tiny_ldt.h"

using namespace std;

Rand jyorandengine;
hitable_list world;

// 拒绝方法随机生成球内一点
vec3 randomInUnitSphere() {
  vec3 p;
  do {
    p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1));
  } while (p.squared_length() >= 1.0);
  return p;
}

// 拒绝方法随机生成半球内一点
vec3 randomInHemisphere(const vec3& normal) {
  // 随机生成一个单位球内的点
  vec3 in_unit_sphere = randomInUnitSphere();
  // 如果点在法线方向的半球上
  if (dot(in_unit_sphere, normal) > 0.0)
    return in_unit_sphere;  // 返回该点作为散射方向
  else
    return -in_unit_sphere;  // 否则返回该点的反方向作为散射方向
}

// 拒绝方法随机生成圆内一点
vec3 randomInUnitDisk() {
  vec3 p;
  do {
    p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
             jyorandengine.jyoRandGetReal<double>(-1, 1), 0);
  } while (p.squared_length() >= 1.0);
  return p;
}

// 反演方法余弦采样生成半球面上一点
vec3 randomCosineDirection() {
  auto r1 = jyorandengine.jyoRandGetReal<double>(0, 1);
  auto r2 = jyorandengine.jyoRandGetReal<double>(0, 1);

  auto phi = 2 * PI * r1;
  auto x = cos(phi) * sqrt(r2);
  auto y = sin(phi) * sqrt(r2);
  auto z = sqrt(1 - r2);

  return vec3(x, y, z);
}

tiny_ldt<float>::light ldt;
vector<vector<float>> intensityDis;

float getIntesiy(float C, float gamma){
  assert(C>=0 && C<=2*M_PI && gamma>=0 && gamma<=M_PI);
  int Cindex = floor(C/M_PI*180.0/ldt.dc);
  int gammaindex = floor(gamma/M_PI*180.0/ldt.dg);
  if(gamma==0.0 || gamma==M_PI)
    return intensityDis[Cindex][gammaindex];


  float d = 0.0, e = 0.0;
  while(d+ldt.dg<=gamma/M_PI*180.0)
    d += ldt.dg;
  while(e+ldt.dc<=C/M_PI*180.0)
    e += ldt.dc;
  float a = 1.0-(gamma/M_PI*180.0-d)/ldt.dg;
  float b = 1.0-(C/M_PI*180.0-e)/ldt.dc;
  float value1 = (a*intensityDis[Cindex][gammaindex]+(1-a)*intensityDis[Cindex][gammaindex+1]);
  float value2 = (a*intensityDis[Cindex+1][gammaindex]+(1-a)*intensityDis[Cindex+1][gammaindex+1]);
  return b*value1 + (1-b)*value2;
}

// 颜色着色
vec3 color(const ray& in, int depth) {
  hit_record rec;
  // 减少误差，-0.00001也可以是交点
  if (world.hitanything(in, 0.001, DBL_MAX, rec)) {
    scatter_record srec;
    vec3 emitted = rec.mat_ptr->emitted(in, rec, rec.u, rec.v, rec.p);
    if (depth < 5 && rec.mat_ptr->scatter(in, rec, srec)) {
      // 金属和玻璃材质跳过pdf
      if (srec.skip_pdf)
        return srec.attenuation * color(srec.out_ray, depth + 1);

      // 散射光线概率密度函数
      double scatteringpdf = rec.mat_ptr->scattering_pdf(in, rec, srec.out_ray);
      // 一些光源照射不到的地方，返回-1
      if (srec.pdf == -1)
        return emitted;
      return emitted + scatteringpdf * srec.attenuation *
                           color(srec.out_ray, depth + 1) / srec.pdf;
    } else {
      // 直视光源则可以看到光源原本的颜色
      // if (!depth) emitted.make_unit_vector();
      vec3 v = unit_vector(-in.direction());
      return emitted*getIntesiy(atan2(-v.y(), -v.z()) + M_PI, M_PI - acos(-v.x()))/abs(dot(unit_vector(-in.direction()), unit_vector(vec3(-1, 0, 0))))/(80 - 30) / (350 - 300);
    }
  } else {
    // 啥也没打到
    return vec3(0, 0, 0);
  }
  exit(0);
}

std::vector<shared_ptr<hitable>> worldlist;
void buildWorld() {
  texture* whitelightptr = new constant_texture(vec3(500, 500, 500));
  texture* mikuptr = new constant_texture(vec3(0.223, 0.773, 0.733));
  texture* redptr = new constant_texture(vec3(0.65, 0.05, 0.05));
  texture* whiteptr = new constant_texture(vec3(0.73, 0.73, 0.73));
  texture* greenptr = new constant_texture(vec3(0.12, 0.45, 0.15));
  texture* groundtexptr = new constant_texture(vec3(0.48, 0.83, 0.53));
  texture* metalptr = new constant_texture(vec3(0.8, 0.85, 0.88));
  texture* checkertextptr =
      new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)),
                          new constant_texture(vec3(0.9, 0.9, 0.9)));

  texture* metaltextureptr = new constant_texture(
      vec3(0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1)),
           0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1)),
           0.5 * (1 + jyorandengine.jyoRandGetReal<double>(0, 1) *
                          jyorandengine.jyoRandGetReal<double>(0, 1))));
  texture* noisetextptr = new noise_texture(0.01);

  int nx, ny, nn;
  unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
  texture* imagetextureptr = new image_texture(tex_data, nx, ny);

  // worldlist.emplace_back(
  //     new rectangle_yz(0, 555, 0, 555, 555, new lambertian(whiteptr)));
  // worldlist.emplace_back(
  //     new rectangle_yz(0, 555, 0, 555, 0, new lambertian(greenptr)));
  // worldlist.emplace_back(new rectangle_xz(213, 343, 227, 332, 554,
  //                                        new diffuse_light(whitelightptr)));
  worldlist.emplace_back(new rectangle_yz(30, 80, 300, 350, 500,
                                         new diffuse_light(whitelightptr)));

  // worldlist.emplace_back(new sphere(vec3(500, 55, 325), 15, new diffuse_light(whitelightptr)));
  // worldlist.emplace_back(
  //     new rectangle_xz(0, 555, 0, 555, 555, new lambertian(whiteptr)));
  worldlist.emplace_back(
     new rectangle_xz(-1000, 1000, -1000, 1000, 0, new lambertian(whiteptr)));
  // worldlist.emplace_back(
  //     new rectangle_xy(0, 555, 0, 555, 555, new lambertian(whiteptr)));

  world = hitable_list(worldlist);
}

int getfileline() {
  std::ifstream file("output.PPM");
  // 判断文件是否打开成功
  if (file.is_open()) {
    int line_count = 0;
    std::string line;
    while (std::getline(file, line)) {
      ++line_count;
    }
    return line_count;
  } else
    return -1;
}

int main() {
  
  std::string err;
  std::string warn;
  if (!tiny_ldt<float>::load_ldt("photometry\\SLOTLIGHT_42184612.LDT", err, warn, ldt)) {
    cout << "failed" << endl;
  }
  if (!err.empty()) 
    cout << err << endl;
  if (!warn.empty())
    cout << warn << endl;
  
  int cnt = 0;
  cout << ldt.dc << endl;//15
  cout << ldt.dg << endl;//5

  for(float i = 0.0; i <= 360.0; i += ldt.dc){
    int j = i;
    // 105/15=7 259/37-1=6
    if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1 &&
    (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc))
    // 105 %= (259/37-1)=6*15
      j %= (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc);
    else if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)
      j = 0;
    if(i==270 && (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc>=90))
      j = 90;
    // cout << "i" << i << " j" << j << endl;
    // cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1) << endl;
    // cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1)+(int)(180.0/ldt.dg)+1 << endl;
    intensityDis.emplace_back(
        vector<float>(ldt.luminous_intensity_distribution.begin()+(int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1)
                    , ldt.luminous_intensity_distribution.begin()+(int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1)+(int)(180.0/ldt.dg)+1)
    );
  }
  for(auto v: intensityDis){
    for(auto p: v)
      cout << p << " ";
    cout << endl;
  }
  cout << intensityDis.size() << endl;
  system("pause");

  // 是否重新渲染
  int startoveragain = 1;

  int curline = getfileline();

  ofstream mout;
  if (startoveragain)
    mout.open("output.PPM");
  else
    mout.open("output.PPM", ios::app);

  // 画布的长
  int nx = 800;
  // 画布的宽
  int ny = 400;
  // 画布某一点的采样数量
  int ns = 100;

  buildWorld();
  // 正常视角
  // vec3 lookfrom(208, 75, -200), lookat(298, 75, 0);
  // 俯视
  vec3 lookfrom(499.99, 800, 325), lookat(500, 0, 324.999999);
  camera cam(lookfrom, lookat, 40, double(nx) / double(ny), 0.0, 10.0, 0.0,
             1.0);

  int pauseflag = 1;
  int si, sj;
  if (curline <= 3 || startoveragain) {
    mout << "P3\n" << nx << " " << ny << "\n255\n";
    sj = ny - 1, si = 0;
  } else if (curline == nx * ny + 3) {
    system("pause");
    return 0;
  } else {
    curline -= 3;
    sj = ny - 1 - curline / nx;
    si = curline - nx * (curline / nx);
  }

  int sqrtns = int(sqrt(ns));
  double resqrtns = 1.0 / sqrtns;
  for (int j = sj; j >= 0; j--) {
    cout << "loading..." << 100 * (1.0 - double(j) / double(ny)) << "%" << endl;
    int starti = pauseflag ? si : 0;
    pauseflag = 0;
    for (int i = starti; i < nx; i++) {
      // 最终的颜色
      vec3 col(0, 0, 0);
      for (int dj = 0; dj < sqrtns; dj++)
        for (int di = 0; di < sqrtns; di++) {
          // 蒙特卡洛-抖动采样，将像素划分成更密的小格子，每个格子里随机取一个点采样
          double uplus =
              -0.5 + resqrtns * ((double)di +
                                 jyorandengine.jyoRandGetReal<double>(-1, 1));
          double vplus =
              -0.5 + resqrtns * ((double)dj +
                                 jyorandengine.jyoRandGetReal<double>(-1, 1));

          // 点(u,v)是点(i,j)的反离散化
          double u = (double(i) + uplus) / double(nx);
          double v = (double(j) + vplus) / double(ny);

          // 一条射向画布上点(u,v)的光线，注意(u,v)不是真实坐标而是在画布上的比例位置
          ray r = cam.get_ray(u, v);
          vec3 cr = color(r, 0);
          cr.e[0] = max(cr.e[0], 0.0), cr.e[0] = min(cr.e[0], 1.0);
          cr.e[1] = max(cr.e[1], 0.0), cr.e[1] = min(cr.e[1], 1.0);
          cr.e[2] = max(cr.e[2], 0.0), cr.e[2] = min(cr.e[2], 1.0);
          col += cr;
        }
      // 取颜色的平均值
      col /= double(ns);
      // 消除坏点
      
      // gamma2修正，提升画面的质量
      col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
      int ir = int(255.99 * col[0]);
      int ig = int(255.99 * col[1]);
      int ib = int(255.99 * col[2]);
      ir = max(0, ir), ig = max(0, ig), ib = max(0, ib);
      ir = min(ir, 255), ig = min(ig, 255), ib = min(ib, 255);
      stringstream ss;
      ss << ir << " " << ig << " " << ib << "\n";
      mout << ss.str();
    }
    system("cls");
  }

  mout.close();
  system("pause");
  return 0;
}