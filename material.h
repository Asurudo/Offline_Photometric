#ifndef MATERIALH
#define MATERIALH

#include <iostream>

#include "hitable.h"
#include "jyorand.h"
#include "onb.h"
#include "pdf.h"
#include "ray.h"
#include "texture.h"

// 前向声明
struct hit_record;
extern vec3 randomInUnitSphere();
extern vec3 randomInHemisphere(const vec3& normal);
extern vec3 randomCosineDirection();
extern Rand jyorandengine;

class scatter_record{
  public: 
    vec3 attenuation;
    double pdf;
    ray out_ray;
    bool skip_pdf;
};

class material {
 public:
  // 材质指针
  texture* textureptr;
  material() {}
  material(texture* a) : textureptr(a) {}
  // 计算反射光线
  virtual ray reflect(const ray& r_in, const hit_record& rec) const = 0;
  // 计算总共的散射光线
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const = 0;
  // 发光函数
  virtual vec3 emitted(const ray& r_in, const hit_record& rec, double u,
                       double v, const vec3& p) const {
    // 返回纯黑，表示不发光
    return vec3(0, 0, 0);
  }

  // 散射概率密度函数
  virtual double scattering_pdf(const ray& r_in, const hit_record& rec,
                                const ray& scattered) const {
    return -1.0;
  }
};

class lambertian : public material {
 public:
  lambertian(texture* a) : material(a) {}
  virtual ray reflect(const ray& r_in, const hit_record& rec) const override {
    vec3 target = rec.p - r_in.A + rec.normal + randomInUnitSphere();
    return ray(rec.p, unit_vector(target - (rec.p - r_in.A)), r_in.time());
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const override {
    // 不考虑pdf的漫反射
    // srec.out_ray = reflect(r_in, rec);
    // srec.attenuation = textureptr->value(rec.u, rec.v, rec.p);
    // srec.skip_pdf = true;
    // return true;
    
    // 光源采样与球体玻璃混合采样
    // mixture_pdf mixed_pdf(std::make_shared<mixture_pdf> (std::make_shared<hitable_pdf>(rec.p, rec.normal),
    //                       std::make_shared<cosine_pdf>(rec.normal), 0.8),
    //                       std::make_shared<sphere_dielectric_pdf>(rec.p, rec.normal), 0.95);
    // srec.out_ray = ray(rec.p, mixed_pdf.generate(), r_in.time());
    // srec.attenuation = textureptr->value(rec.u, rec.v, rec.p);
    // srec.pdf = mixed_pdf.value(srec.out_ray.direction());
    // srec.skip_pdf = false;
    sphere_dielectric_pdf sdp(rec.p, rec.normal);
    // hitable_pdf light_pdf(rec.p, rec.normal);
    srec.out_ray = ray(rec.p, sdp.generate(), r_in.time());
    srec.attenuation = textureptr->value(rec.u, rec.v, rec.p);
    srec.pdf = sdp.value(srec.out_ray.direction());
    srec.skip_pdf = false;
    
    // 球体玻璃采样
    // sphere_dielectric_pdf sdpdf(rec.p, rec.normal);
    // srec.out_ray = ray(rec.p, sdpdf.generate(), r_in.time());
    // srec.attenuation = textureptr->value(rec.u, rec.v, rec.p);
    // srec.pdf = sdpdf.value(srec.out_ray.direction());
    // srec.skip_pdf = false;
    // // 射不到球体就别射了
    // if(srec.pdf <= 0.0)
    //     srec.out_ray = reflect(r_in, rec), srec.skip_pdf = true;
    
    return true;
  }

  virtual double scattering_pdf(const ray& r_in, const hit_record& rec,
                                const ray& scattered) const override {
    // 表面法线与散射光线夹角的余弦值
    double cos_theta =
        dot(unit_vector(rec.normal), unit_vector(scattered.direction()));
    // 如果小于零，代表光线往物体内部射，因此设为0
    return (cos_theta < 0 ? 0 : cos_theta) / PI;
  }
};

class metal : public material {
 public:
  // 模糊镜面反射特有的模糊系数
  double fuzz;
  // 对于rgb各个分量的反射量以及模糊镜面反射的系数[0,1]
  metal(texture* a, double f = 0.75) : material(a) {
    if (f < 1)
      fuzz = f;
    else
      fuzz = 1;
  }
  virtual ray reflect(const ray& r_in, const hit_record& rec) const override {
    return ray(rec.p,
               unit_vector(unit_vector(rec.p - r_in.A) +
                           (2 * dot(-unit_vector(rec.p - r_in.A), rec.normal) *
                            rec.normal) +
                           fuzz * randomInUnitSphere()),
               r_in.time());
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const override {
    srec.out_ray = reflect(r_in, rec);
    srec.attenuation = textureptr->value(0, 0, rec.p);
    srec.skip_pdf = true;
    return true;
  }
};

class dielectric : public material {
 public:
  // 折射率
  double refIdx;

  dielectric(double ri) : refIdx(ri) {}

  // 一种简单的方法确定光的反射量，返回一个[0,1]的数值，将其看作一个概率
  double schlick(double cosine, double RI) const {
    double R0 = (1.0 - RI) / (1.0 + RI);
    R0 = R0 * R0;
    return R0 + (1 - R0) * pow((1 - cosine), 5);
  }
  // 镜面反射
  virtual ray reflect(const ray& r_in, const hit_record& rec) const override {
    return ray(rec.p,
               unit_vector(unit_vector(rec.p - r_in.A) +
                           (2 * dot(-unit_vector(rec.p - r_in.A), rec.normal) *
                            rec.normal)),
               r_in.time());
  }
  bool refract(const ray& r_in, const hit_record& rec, const vec3& n,
               double niOverNt, ray& refracted) const {
    vec3 uIn = unit_vector(rec.p - r_in.A);
    double dt = dot(uIn, n);
    double delta = 1.0 - niOverNt * niOverNt * (1 - dt * dt);

    bool returnValue;
    refracted = ((returnValue = (delta > 0))
                     ? ray(rec.p, niOverNt * (uIn - n * dt) - n * sqrt(delta),
                           r_in.time())
                     : ray());

    return returnValue;
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const override {
    // 入射介质的折射率比上折射介质的折射率
    double niOverNt, cosine, reflectProb = 1.0;
    // 折射介质的法线
    vec3 outwardNormal;
    // 折射光线与反射光线
    ray refracted, reflected = reflect(r_in, rec);
    // 反射率永远为1，因为折射的电解质材料不吸收任何光线，全部反射出去
    srec.attenuation = vec3(1.0, 1.0, 1.0);

    // 光线和法线呈锐角
    if (dot(r_in.direction(), rec.normal) > 0) {
      outwardNormal = -rec.normal;
      niOverNt = refIdx / 1.0;
      cosine = refIdx * dot(r_in.direction(), rec.normal) /
               r_in.direction().length();
    } else {
      outwardNormal = rec.normal;
      niOverNt = 1.0 / refIdx;
      cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
    }
    // 如果折射(未发生全反射现象，也就是入射角小于临界角)，返回折射光线，否则返回反射光线
    if (refract(r_in, rec, outwardNormal, niOverNt, refracted))
      reflectProb = schlick(cosine, refIdx);

    // 将反射量看作概率，随机一下，如果落在概率内，则返回反射光线，否则返回折射光线
    srec.out_ray = ((jyorandengine.jyoRandGetReal<double>(0, 1) < reflectProb)
                     ? reflected
                     : refracted);
    srec.skip_pdf = true;
    return true;
  }
};

class diffuse_light : public material {
 public:
  diffuse_light(texture* a) : material(a) {}
  virtual ray reflect(const ray& r_in, const hit_record& rec) const override {
    return ray(vec3(0, 0, 0), vec3(39, 39, 39));
    exit(0);
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const override {
    return false;
  }
  virtual vec3 emitted(const ray& r_in, const hit_record& rec, double u,
                       double v, const vec3& p) const override {
    // 保证光源只向下发射光线
    // if (dot(vec3(0, -1, 0), r_in.direction()) < 0.0)
      return textureptr->value(u, v, p);
    // return vec3(0, 0, 0);
  }
};

// 各向同性
class isotropic : public material {
 public:
  isotropic(texture* a) : material(a) {}
  virtual ray reflect(const ray& r_in, const hit_record& rec) const override {
    return ray(vec3(0, 0, 0), vec3(39, 39, 39));
    exit(0);
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       scatter_record& srec) const override {
    // 往任意方向反射光线
    // 和粗糙磨砂表面的区别是，粗糙磨砂表面不会往物体内反射
    srec.out_ray = ray(rec.p, randomInUnitSphere());
    srec.attenuation = textureptr->value(rec.u, rec.v, rec.p);
    srec.pdf = 1 / (4 * PI);
    srec.skip_pdf = false;
    return true;
  }

  virtual double scattering_pdf(const ray& r_in, const hit_record& rec,
                                const ray& scattered) const override {
    return 1 / (4 * PI);
  }
};

#endif