#ifndef MATERIALH
#define MATERIALH

#include <iostream>

#include "hitable.h"
#include "jyorand.h"
#include "ray.h"
#include "texture.h"

// 前向声明
struct hit_record;
extern vec3 randomInUnitSphere();
extern Rand jyorandengine;

class material {
 public:
  // 材质指针
  texture* textureptr;
  material() {}
  material(texture* a) : textureptr(a) {}
  // 计算反射光线
  virtual ray reflect(const ray& r_in, const hit_record& rec) const = 0;
  virtual bool scatter(const ray& r_in, const hit_record& rec,
                       vec3& attenuation, ray& scattered) const = 0;
  virtual vec3 emitted(double u, double v, const vec3& p) const {
    // 返回纯黑，表示不发光
    return vec3(0, 0, 0);
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
                       vec3& attenuation, ray& scattered) const override {
    scattered = reflect(r_in, rec);
    attenuation = textureptr->value(rec.u, rec.v, rec.p);
    return true;
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
                       vec3& attenuation, ray& scattered) const override {
    scattered = reflect(r_in, rec);
    attenuation = textureptr->value(0, 0, rec.p);
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
                       vec3& attenuation, ray& scattered) const override {
    // 入射介质的折射率比上折射介质的折射率
    double niOverNt, cosine, reflectProb = 1.0;
    // 折射介质的法线
    vec3 outwardNormal;
    // 折射光线与反射光线
    ray refracted, reflected = reflect(r_in, rec);
    // 反射率永远为1，因为折射的电解质材料不吸收任何光线，全部反射出去
    attenuation = vec3(1.0, 1.0, 1.0);

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
    scattered = ((jyorandengine.jyoRandGetReal<double>(0, 1) < reflectProb)
                     ? reflected
                     : refracted);

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
                       vec3& attenuation, ray& scattered) const override {
    return false;
  }
  virtual vec3 emitted(double u, double v, const vec3& p) const override {
    // 返回材质颜色
    return textureptr->value(u, v, p);
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
                       vec3& attenuation, ray& scattered) const override {
    // 往任意方向反射光线
    // 和粗糙磨砂表面的区别是，粗糙磨砂表面不会往物体内反射
    scattered = ray(rec.p, randomInUnitSphere());
    attenuation = textureptr->value(rec.u, rec.v, rec.p);
    return true;
  }
};

#endif