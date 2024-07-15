#ifndef PDFH
#define PDFH

#include "hitable.h"
#include "onb.h"

extern vec3 randomInUnitSphere();
extern vec3 randomInHemisphere(const vec3& normal);
extern vec3 randomCosineDirection();
extern Rand jyorandengine;

class pdf {
 public:
  // 分母上的pdf，也就是采样服从的pdf的值，direction是generate的光线方向
  virtual double value(const vec3& direction) = 0;
  // 按照pdf采样出来的光线
  virtual vec3 generate() = 0;
};

// 均匀采样-半球表面散射
class sphere_pdf : public pdf {
 private:
  vec3 normal;

 public:
  sphere_pdf(vec3 normal) : normal(normal) {}

  virtual double value(const vec3& direction) override { return 1 / (4 * PI); }

  virtual vec3 generate() override { return randomInHemisphere(normal); }
};

// 余弦采样-半球表面散射
class cosine_pdf : public pdf {
 private:
  // 以法线为一轴创建一个坐标系
  onb uvw;

 public:
  // 传入法线
  cosine_pdf(const vec3& w) : uvw(w) {}

  virtual double value(const vec3& direction) override {
    // 单位半球面上按照余弦采样生成的点转化为相对坐标
    double cosine_theta = dot(unit_vector(uvw.w()), unit_vector(direction));
    return std::max(0.0, cosine_theta / PI);
  }

  virtual vec3 generate() override {
    vec3 scatter_direction = uvw.local(randomCosineDirection());

    if (scatter_direction.near_zero()) scatter_direction = uvw.w();

    return unit_vector(scatter_direction);
  }
};

// 光源采样
class hitable_pdf : public pdf {
 private:
  // 物体与光线的交点
  vec3 p;
  // 接触点到光源上一点的平方
  double distance_squard;
  // 法线
  vec3 normal;

 public:
  hitable_pdf(const vec3& p, const vec3& normal) : p(p), normal(normal) {}

  virtual double value(const vec3& direction) override {
    // 夹角为钝角，光源在物体背面，届不到
    if (dot(direction, normal) < 0) return -1;
    // 光源面积
    double light_area = (85 - 30) * (355 - 300);
    // 单位向量巧妙计算夹角
    double light_cosine = fabs(direction.y());
    // 如果光源和物体表面平行，也算届不到
    if (light_cosine < 1e-5) return -1;

    return distance_squard / (light_area * light_cosine);
  }

  virtual vec3 generate() override {
    // 强势硬编码，接管光线处理
    // 光源上随机一点
    vec3 on_light = vec3(500, jyorandengine.jyoRandGetReal<double>(30, 85),
                         jyorandengine.jyoRandGetReal<double>(300, 355));
    // 接触点到光源的射线
    vec3 to_light = on_light - p;
    // 计算接触点到光源上一点的平方
    distance_squard = to_light.squared_length();
    return to_light = unit_vector(to_light);
  }
};

class sphere_dielectric_pdf : public pdf{
  private:
    vec3 p;
    vec3 normal;
    // 球体硬编码
    vec3 center;
    double radius;
    // 均匀采样从任意点可见的球体的一侧
    vec3 random_to_sphere(double distace_squared){
      double r1 = jyorandengine.jyoRandGetReal<double>(0, 1);
      double r2 = jyorandengine.jyoRandGetReal<double>(0, 1);
      double z = 1 + r2 * (sqrt(1-radius*radius/distace_squared) - 1);

      double phi = 2*PI*r1;
      double x = cos(phi)*sqrt(1-z*z);
      double y = sin(phi)*sqrt(1-z*z);

      return vec3(x, y, z);
    }
  public:
   sphere_dielectric_pdf(const vec3& p, const vec3& normal) : normal(normal), p(p), center(500, 55, 325), radius(15.0) {}

  virtual double value(const vec3& direction) override {
    // 夹角为钝角，球体在物体背面，届不到
    if (dot(direction, normal) < 0) return -1;
    
    double cos_theta_max = sqrt(1 - radius*radius/(center-p).squared_length());
    double solid_angle = 2*PI*(1-cos_theta_max);

    return 1.0 / solid_angle;
  }

  virtual vec3 generate() override {
    vec3 direction = center - p;
    double distance_squared = direction.squared_length();
    onb uvw(direction);
    return unit_vector(uvw.local(random_to_sphere(distance_squared)));
  }
};

class mixture_pdf : public pdf {
 private:
  std::shared_ptr<pdf> p[2];
  double zeropro;
  int tp;

 public:
  mixture_pdf(std::shared_ptr<pdf> p0, std::shared_ptr<pdf> p1,
              double pro = 0.5) {
    p[0] = p0;
    p[1] = p1;
    zeropro = pro;
  }

  virtual double value(const vec3& direction) override {
    return p[tp]->value(direction);
  }

  virtual vec3 generate() override {
    tp = jyorandengine.jyoRandGetReal<double>(0, 1) > zeropro;
    return p[tp]->generate();
  }
};

#endif