//#define COSINE_SAMPLING
//#define LIGHT_SAMPLING
#define LIGHT_DOUBLEAXIS_SAMPLE

// 画布的长
int nx = 800;
// 画布的宽
int ny = 600;
// 画布某一点的采样数量
int ns = 1500;


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
#include "perlin.h"
#include "rectangle.h"
#include "smoke.h"
#include "sphere.h"
#include "stb_image.h"
#include "texture.h"
#include "transformation.h"
#include "tiny_ldt.h"

using namespace std;

// double axis moment
int n1 = 5;
int n2 = 1;
vec3 axis_w(0,1,0); 
vec3 axis_v(0,1,0);

std::string filename = "PANOS_60813872.LDT";
double roughness = 0.05;
//vec3 lookfrom(0, 60, 0), lookat(0.0001, 0, 0);
// vec3 lookfrom(25, 15, 20), lookat(0, 0, 0.029);
//vec3 lookfrom(25, 2, 0), lookat(0, 2, 0);
vec3 lookfrom(-10, 3, 0), lookat(5, 1, 0);

Rand jyorandengine;
hitable_list world;

vec3 randomInUnitSphere() {
  // vec3 p;
  // do {
  //   p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
  //            jyorandengine.jyoRandGetReal<double>(-1, 1),
  //            jyorandengine.jyoRandGetReal<double>(-1, 1));
  // } while (p.squared_length() >= 1.0);
  // return p;

    // const double xi1 = jyorandengine.jyoRandGetReal<double>(0, 1);
		// const double xi2 = jyorandengine.jyoRandGetReal<double>(0, 1);
		// const double theta_h = acos( sqrt(1-xi1) );
		// const double cosph = cos( 2.0 * M_PI * xi2 );
		// const double sinph = sin( 2.0 * M_PI * xi2 );
		// const double costh = cos( theta_h );
		// const double sinth = sin( theta_h );
    // const double checkVal = sinth*cosph*sinth*cosph +  sinth*sinph*sinth*sinph +  costh*costh;
    // assert(checkVal >=0.99  && checkVal <= 1.01);
    // return vec3( sinth * cosph, costh, sinth * sinph );

    vec3 p = randomInUnitDisk();
    p.e[1] = sqrt(1.0 - p.x()*p.x() - p.z()*p.z());
    return p;
}

vec3 randomInUnitDisk() {
  vec3 p;
  do {
    p = vec3(jyorandengine.jyoRandGetReal<double>(-1, 1),
             0, jyorandengine.jyoRandGetReal<double>(-1, 1));
  } while (p.squared_length() >= 1.0);
  return p;
}

tiny_ldt<float>::light ldt;
vector<vector<float>> intensityDis;

float getIntesiy(float C, float gamma){
  assert(C>=0 && C<=2*M_PI && gamma>=0 && gamma<=M_PI);
  int Cindex = floor(C/M_PI*180.0/ldt.dc);
  int gammaindex = floor(gamma/M_PI*180.0/ldt.dg);
  //if(gamma==0.0 || gamma==M_PI)
  //  return intensityDis[Cindex][gammaindex];

  float d = 0.0, e = 0.0;
  while(d+ldt.dg<=gamma/M_PI*180.0)
    d += ldt.dg;
  while(e+ldt.dc<=C/M_PI*180.0)
    e += ldt.dc;
  float a = 1.0-(gamma/M_PI*180.0-d)/ldt.dg;
  float b = 1.0-(C/M_PI*180.0-e)/ldt.dc;
  float value1 = (a*intensityDis[Cindex][gammaindex]+(1-a)*intensityDis[Cindex][gammaindex+1]);
  float value2 = (a*intensityDis[Cindex+1][gammaindex]+(1-a)*intensityDis[Cindex+1][gammaindex+1]);
  return 600 * (b*value1 + (1-b)*value2)/683.f*10;
}

vec3 m_t, m_b, m_n;

//(0,1,0)
void getMTNB(const vec3 &n)
{
  m_n = n; 
  m_t = unit_vector(cross(n, (std::abs(n.x()) < std::abs(n.z())) ? vec3(0.0, -n.z(), n.x()) : vec3(-n.y(), n.x(), 0.0)));
  m_b = cross(m_t, m_n);
}

vec3 to_local(const vec3 &w)
{
  // m_n = n = (0, 1, 0)
  // return vec3(dot(m_t, w), dot(m_b, w), dot(m_n, w)); 
  return vec3(dot(m_t, w), dot(m_n, w), dot(m_b, w)); 
}

double d(const vec3 &lwi, const vec3 &lwo, const double &m_alpha)
{
	const vec3 h = unit_vector( lwi + lwo );
  const auto cos2 = h.y() * h.y();
	const auto sin2 = max(0.0, 1.0 - cos2);
  
  const double denom = M_PI * m_alpha * m_alpha * ( cos2 + sin2 / ( m_alpha * m_alpha ) ) * ( cos2 + sin2 / ( m_alpha * m_alpha ) );
  assert( std::isfinite( denom ) );
  return 1.0 / denom;

  // return m_alpha * m_alpha / ( M_PI * ( cos2 * m_alpha * m_alpha + sin2 ) * ( cos2 * m_alpha * m_alpha + sin2 ) );
}

//幾何項を返す
double g(const vec3 &lwi, const vec3 &lwo, const double &m_alpha)
{
	const double tan_i = 1.0 / ( lwi.y() * lwi.y() ) - 1.0; assert( std::isfinite( tan_i ) );
	const double tan_o = 1.0 / ( lwo.y() * lwo.y() ) - 1.0; assert( std::isfinite( tan_o ) );
  
	const double lambda_i = ( - 1.0 + sqrt( 1.0 + m_alpha * m_alpha * tan_i ) ) / 2.0; assert( std::isfinite( lambda_i ) );
	const double lambda_o = ( - 1.0 + sqrt( 1.0 + m_alpha * m_alpha * tan_o ) ) / 2.0; assert( std::isfinite( lambda_o ) );
	
  return 1.0 / ( 1.0 + lambda_i + lambda_o );
}

double f(const vec3& lwi, const vec3& lwo, double f0)
{
    const vec3 h = unit_vector( lwi + lwo );
    const double cosine = dot( h, lwi );
    const double tmp = ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine ) * ( 1.0 - cosine );
    return f0 + ( 1.0 - f0 ) * tmp;
}


float eval(const vec3& V, const vec3& L, const float alpha)
{
		if(V.y() <= 0)
		{
      std::cout << "V.y() <= 0" << std::endl;
      return 0;
		}

		// masking
		const float a_V = 1.0f / alpha / tanf(acosf(V.y()));
		const float LambdaV = (V.y()<1.0f) ? 0.5f * (-1.0f + sqrtf(1.0f + 1.0f/a_V/a_V)) : 0.0f;
    //const float LambdaV = (V.y()<1.0f) ? (1.0f - 1.259f*a_V + 0.396f*a_V*a_V) / (3.535f*a_V + 2.181f*a_V*a_V) : 0.0f;
		const float G1 = 1.0f / (1.0f + LambdaV);

		// shadowing
		float G2;
		if(L.y() <= 0.0f)
			G2 = 0;
		else
		{
			const float a_L = 1.0f / alpha / tanf(acosf(L.y()));
			const float LambdaL = (L.y()<1.0f) ? 0.5f * (-1.0f + sqrtf(1.0f + 1.0f/a_L/a_L)) : 0.0f;
      //const float LambdaL = (L.y()<1.0f) ? (1.0f - 1.259f*a_L + 0.396f*a_L*a_L) / (3.535f*a_L + 2.181f*a_L*a_L) : 0.0f;
			G2 = 1.0f / (1.0f + LambdaV + LambdaL);
		}

		// D
		const vec3 H = unit_vector(V+L);
		const float slopex = H.x()/H.y();
		const float slopez = H.z()/H.y();
		float D = 1.0f / (1.0f + (slopex*slopex+slopez*slopez)/alpha/alpha);
		D = D*D;
		D = D / (3.14159f * alpha * alpha * H.y()*H.y()*H.y()*H.y());
    //float D = expf(-(slopex*slopex + slopez*slopez)/(alpha*alpha)) / (3.14159f * alpha * alpha * H.y()*H.y()*H.y()*H.y());

		float res = D * G2 / 4.0f / V.y();
    //cout << res << endl;
		return res;
}

// マクロ法線 光源L 照相机V ラフネス値 f0
double BRDF_Specular_GGX(vec3 N, vec3 L, vec3 V, double roughness, double f0)
{
    L = to_local(L);
    V = to_local(V);
    //std::cout << eval(L, V, roughness) << std::endl;
    return eval(V, L, roughness);
    // vec3 H = unit_vector(L + V);
    double NoV = V.y();// dot(N, V);
    double NoL = L.y();// dot(N, L);
    
    // NoV = abs(NoV), NoL = abs(NoL);

    // if(NoV <= 0.0 || NoL <= 0.0)
    //  return 0.0;

    // double NoH = dot(N, H);
    // double LoH = dot(L, H);
    
    // 法線分布
    double D = d(L, V, roughness);
    // 幾何減衰
    double G = g(L, V, roughness);
    // フレネル
    double F = f(L, V, f0);
    // cout << D * G * F / (4.0 * NoL * NoV) << endl;
    // スペキュラBRDF
    return D * G * F / (4.0 * NoL * NoV);
}

// 颜色着色
vec3 color(const ray& in, int depth) {
  
  #ifdef LIGHT_SAMPLING
    assert(depth<2);
  #endif
  
  hit_record rec;
  if (world.hitanything(in, 0.0001, DBL_MAX, rec)) {
    
    // if(rec.p.x()<=-0.0001)
    //   return vec3(0, 0, 0);
    
    // 反射光
    ray scattered;
    // 吸收度
    vec3 attenuation;
    vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    if (depth < 5 && rec.mat_ptr->scatter(in, rec, attenuation, scattered)){
      // 余弦
      double cos_theta = dot(unit_vector(rec.normal), unit_vector(scattered.direction()));
      // double brdf = 1.0 / PI;
      double brdf = (n1+2)/(2*PI);
      assert(rec.normal.x()==0 && rec.normal.y()==1 && rec.normal.z()==0);
      // double brdf = BRDF_Specular_GGX(unit_vector(rec.normal), 
      //                                 unit_vector(scattered.direction()), 
      //                                 unit_vector(-in.direction()), 
      //                                 roughness, 1.0); 

      #ifdef COSINE_SAMPLING
      return brdf * PI * color(scattered, depth + 1);
      #endif

      #ifdef LIGHT_SAMPLING
      return brdf * max(cos_theta, 0.0)* color(scattered, depth + 1);
      #endif

      #ifdef LIGHT_DOUBLEAXIS_SAMPLE
      return brdf * color(scattered, depth + 1);
      #endif
    }
    else {
      // 光源
      if (!depth) return vec3(1, 1, 1);
      vec3 v = unit_vector(-in.direction());
      #ifdef COSINE_SAMPLING
      if(dot(unit_vector(-in.direction()), vec3(1, 0, 0))>0)
      return emitted*getIntesiy(atan2(-v.y(), -v.z()) + M_PI, M_PI - acos(-v.x()))
        /dot(unit_vector(-in.direction()), vec3(1, 0, 0)) / (3.4-0.4) / (3.0-0.0);
      else
        return emitted*getIntesiy(atan2(-v.y(), -v.z()) + M_PI, M_PI - acos(-v.x()))
        /dot(unit_vector(-in.direction()), vec3(-1, 0, 0)) / (3.4-0.4) / (3.0-0.0);
      #endif

      #ifdef LIGHT_SAMPLING
      return emitted*getIntesiy(atan2(-v.y(), -v.z()) + M_PI, M_PI - acos(-v.x()))
                                /dot(rec.p-in.origin(), rec.p-in.origin()); // I/distance^2
      #endif

      #ifdef LIGHT_DOUBLEAXIS_SAMPLE
      vec3 p2q = rec.p - in.origin();
      axis_w = unit_vector(axis_w);
      axis_v = unit_vector(axis_v);
      double dam = pow(dot(axis_w, unit_vector(p2q)), n1) * pow(dot(axis_v, unit_vector(p2q)), n2);
      double cos_theta_prime = dot(-unit_vector(p2q), vec3(0, -1, 0));
      if(cos_theta_prime < 0.0)
        cos_theta_prime = 0.0;
      if(dam < 0.0)
        dam = 0.0;
      assert(cos_theta_prime >= 0.0 && dam >= 0.0);
      double rnt = ((1.5-(-1.5)) * (1.5-(-1.5))*dam*cos_theta_prime)/(dot(rec.p-in.origin(), rec.p-in.origin()));
      return vec3(rnt, rnt, rnt);
      #endif
    }
  } else {
    return vec3(0, 0, 0); // 闇に射る
  }
  exit(0);
}
std::vector<shared_ptr<hitable>> worldlist;
void buildWorld() {
  getMTNB(vec3(0,1,0));
  texture* whitelightptr = new constant_texture(vec3(1, 1, 1));
  texture* mikulightptr = new constant_texture(vec3(0.223, 0.773, 0.733) * 15);
  texture* mikuptr = new constant_texture(vec3(0.223, 0.773, 0.733));
  texture* redptr = new constant_texture(vec3(0.65, 0.05, 0.05));
  texture* whiteptr = new constant_texture(vec3(1, 1, 1));
  texture* greenptr = new constant_texture(vec3(0.12, 0.45, 0.15));
  texture* groundtexptr = new constant_texture(vec3(0.48, 0.83, 0.53));

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

  // // 灯
  // worldlist.emplace_back(new rectangle_yz(0.4, 3.4, -1.5, 1.5, 0,
  //                                        new diffuse_light(whiteptr)));

  worldlist.emplace_back(
    new rectangle_xz(-1.5, 1.5, -1.5, 1.5, 1.0, new diffuse_light(whiteptr))
  );
                                       
  worldlist.emplace_back(
     new rectangle_xz(-40, 40, -40, 40, 0, new lambertian(whiteptr)));

  // 一个玻璃球与一团玻璃球形状的烟雾
  // hitable* glasssphereptr =
  //     new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));
  // worldlist.emplace_back(glasssphereptr);
  // worldlist.emplace_back(new smoke(glasssphereptr, 0.2,
  //                                  new constant_texture(vec3(0.2, 0.4, 0.9))));

  world = hitable_list(worldlist);
}

int getfileline(string filename) {
  std::ifstream file(filename + ".PPM");
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

  if (!tiny_ldt<float>::load_ldt("photometry\\" + filename, err, warn, ldt)) {
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
    // int j = i;
    // // 105/15=7 259/37-1=6
    // if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1 &&
    // (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc))
    // // 105 %= (259/37-1)=6*15
    //   j %= (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc);
    // else if((i/ldt.dc)>ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)
    //   j = 0;
    // if(i==270 && (int)((ldt.luminous_intensity_distribution.size()/((int)(180.0/ldt.dg)+1)-1)*ldt.dc>=90))
    //   j = 90;
    int sz = (180/ldt.dg)+1;
    // cout << "i" << i << " j" << j << endl;
    // cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1) << endl;
    // cout << (int)(j/ldt.dc)*((int)(180.0/ldt.dg)+1)+(int)(180.0/ldt.dg)+1 << endl;
    int st = (int)(i/ldt.dc);
    if(st*sz >= ldt.luminous_intensity_distribution.size())
      st %= ldt.luminous_intensity_distribution.size()/sz;
    intensityDis.emplace_back(
        vector<float>(ldt.luminous_intensity_distribution.begin()+ st*(sz)
                    , ldt.luminous_intensity_distribution.begin()+ st*(sz)+sz)
    );
  }
  for(auto& v: intensityDis){
    for(auto& p: v)
       //if(jyorandengine.jyoRandGetBool(0.01))
        // p += 100;
      cout << p << " ";
    cout << endl;
  }
  cout << intensityDis.size() << endl;
  system("pause");
  
  
  // 是否重新渲染
  int startoveragain = 1;

  int curline = getfileline(filename);

  ofstream mout;
  filename += to_string(roughness);
  if (startoveragain)
    mout.open(filename + ".PPM");
  else
    mout.open(filename + ".PPM", ios::app);


  buildWorld();
  
  camera cam(lookfrom, lookat, 45, double(nx) / double(ny), 0.0, 0.0);

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

  //-----------------------------------------------------------------------
  int sqrtns = int(sqrt(ns));
  double resqrtns = 1.0 / sqrtns;
  for (int j = sj; j >= 0; j--) {
    cout << "loading..." << 100 * (1.0 - double(j) / double(ny)) << "%";
    int starti = pauseflag ? si : 0;
    pauseflag = 0;
    for (int i = starti; i < nx; i++) {
      // 最终的颜色
      vec3 col(0, 0, 0);
      for (int dj = 0; dj < sqrtns; dj++)
        for (int di = 0; di < sqrtns; di++) {
          // 蒙特卡洛-抖动采样，将像素划分成更密的小格子，每个格子里随机取一个点采样
          // double uplus = -0.5 + resqrtns * ((double)di + jyorandengine.jyoRandGetReal<double>(-1, 1)); 
          // double vplus = -0.5 + resqrtns * ((double)dj + jyorandengine.jyoRandGetReal<double>(-1, 1));

          // 点(u,v)是点(i,j)的反离散化
          double u = double(i) / double(nx);
          double v = double(j) / double(ny);

          // 一条射向画布上点(u,v)的光线，注意(u,v)不是真实坐标而是在画布上的比例位置
          ray r = cam.get_ray(u, v);
          col += color(r, 0);
        }
      // 取颜色的平均值
      col /= double(ns);

      // gamma修正，提升画面的质量
      col = vec3(pow(col[0], 1.0/2.2), pow(col[1], 1.0/2.2), pow(col[2], 1.0/2.2));
      int ir = int(255.99 * col[0]);
      int ig = int(255.99 * col[1]);
      int ib = int(255.99 * col[2]);
      ir = min(ir, 255), ig = min(ig, 255), ib = min(ib, 255);
      assert(ir>=0 && ig>=0 && ib>=0);
      // ir = max(ir, 0), ig = max(ig, 0), ib = max(ib, 0);
      stringstream ss;
      ss << ir << " " << ig << " " << ib << "\n";
      mout << ss.str();
    }
    system("cls");
  }

  //-----------------------------------------------------------------------

  // for (int j = sj; j >= 0; j--) {
  //   cout << "loading..." << 100 * (1.0 - double(j) / double(ny)) << "%";
  //   int starti = pauseflag ? si : 0;
  //   pauseflag = 0;
  //   for (int i = starti; i < nx; i++) {
  //     // 最终的颜色
  //     vec3 col(0, 0, 0);
  //     for (int k = 0; k < ns; k++) {
  //       // 点(u,v)是点(i,j)的反离散化
  //       double u = double(i + jyorandengine.jyoRandGetReal<double>(-1, 1)) /
  //                  double(nx);
  //       double v = double(j + jyorandengine.jyoRandGetReal<double>(-1, 1)) /
  //                  double(ny);

  //       //
  //       一条射向画布上点(u,v)的光线，注意(u,v)不是真实坐标而是在画布上的比例位置
  //       ray r = cam.get_ray(u, v);
  //       col += color(r, 0);
  //     }
  //     // 取颜色的平均值
  //     col /= double(ns);
  //     // gamma2修正，提升画面的质量
  //     col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
  //     int ir = int(255.99 * col[0]);
  //     int ig = int(255.99 * col[1]);
  //     int ib = int(255.99 * col[2]);
  //     ir = min(ir, 255), ig = min(ig, 255), ib = min(ib, 255);
  //     stringstream ss;
  //     ss << ir << " " << ig << " " << ib << "\n";
  //     mout << ss.str();
  //   }
  //   system("cls");
  // }

  mout.close();
  system("pause");
  return 0;
}