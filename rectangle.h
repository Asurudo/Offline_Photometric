#ifndef RECTANGLEH
#define RECTANGLEH

#include "material.h"

class rectangle : public hitable {
 protected:
  double dimention_one0, dimention_one1, dimention_two0, dimention_two1,
      dimention_three;
  material* matptr;

 public:
  rectangle() {}
  rectangle(double dimention_one0, double dimention_one1,
                  double dimention_two0, double dimention_two1,
                  double dimention_three, material* matptr)
      : dimention_one0(dimention_one0),
        dimention_one1(dimention_one1),
        dimention_two0(dimention_two0),
        dimention_two1(dimention_two1),
        dimention_three(dimention_three),
        matptr(matptr){};

  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const;
  virtual std::shared_ptr<aabb> getaabb() const;
};

class rectangle_xy : public rectangle {
 public:
  rectangle_xy() {}
  rectangle_xy(double x0, double x1, double y0, double y1, double z,
                     material* matptr)
      : rectangle(x0, x1, y0, y1, z, matptr){};
  std::shared_ptr<aabb> getaabb() const override {
    return std::make_shared<aabb>(
        vec3(dimention_one0, dimention_two0, dimention_three - 0.0001),
        vec3(dimention_one1, dimention_two1, dimention_three + 0.0001));
  };
  bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
    // 联立公式计算出射线的t
    double t = (dimention_three - r.origin().z()) / r.direction().z();

    if (t < tmin || t > tmax) return false;

    // 计算具体的击中坐标
    double x = r.origin().x() + t * r.direction().x();
    double y = r.origin().y() + t * r.direction().y();
    if (x < dimention_one0 || x > dimention_one1 || y < dimention_two0 ||
        y > dimention_two1)
      return false;

    rec.u = (x - dimention_one0) / (dimention_one1 - dimention_one0);
    rec.v = (y - dimention_two0) / (dimention_two1 - dimention_two0);
    rec.t = t;
    rec.mat_ptr = matptr;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 0, 1);
    if (dot(r.direction(), rec.normal) > 0) rec.normal = -rec.normal;
    return true;
  }
};

class rectangle_xz : public rectangle {
 public:
  rectangle_xz() {}
  rectangle_xz(double x0, double x1, double z0, double z1, double y,
                     material* matptr)
      : rectangle(x0, x1, z0, z1, y, matptr){};
  std::shared_ptr<aabb> getaabb() const override {
    return std::make_shared<aabb>(
        vec3(dimention_one0, dimention_three - 0.0001, dimention_two0),
        vec3(dimention_one1, dimention_three + 0.0001, dimention_two1));
  };
  bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
    // 联立公式计算出射线的t
    double t = (dimention_three - r.origin().y()) / r.direction().y();

    if (t < tmin || t > tmax) return false;

    // 计算具体的击中坐标
    double x = r.origin().x() + t * r.direction().x();
    double z = r.origin().z() + t * r.direction().z();
    if (x < dimention_one0 || x > dimention_one1 || z < dimention_two0 ||
        z > dimention_two1)
      return false;

    rec.u = (x - dimention_one0) / (dimention_one1 - dimention_one0);
    rec.v = (z - dimention_two0) / (dimention_two1 - dimention_two0);
    rec.t = t;
    rec.mat_ptr = matptr;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 1, 0);
    if (dot(r.direction(), rec.normal) > 0) rec.normal = -rec.normal;
    return true;
  }
};

class rectangle_yz : public rectangle {
 public:
  rectangle_yz() {}
  rectangle_yz(double y0, double y1, double z0, double z1, double x,
                     material* matptr)
      : rectangle(y0, y1, z0, z1, x, matptr){};
  std::shared_ptr<aabb> getaabb() const override {
    return std::make_shared<aabb>(
        vec3(dimention_three - 0.0001, dimention_one0, dimention_two0),
        vec3(dimention_three + 0.0001, dimention_one1, dimention_two1));
  };
  bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
    // 联立公式计算出射线的t
    double t = (dimention_three - r.origin().x()) / r.direction().x();

    if (t < tmin || t > tmax) return false;

    // 计算具体的击中坐标
    double y = r.origin().y() + t * r.direction().y();
    double z = r.origin().z() + t * r.direction().z();
    if (y < dimention_one0 || y > dimention_one1 || z < dimention_two0 ||
        z > dimention_two1)
      return false;

    rec.u = (y - dimention_one0) / (dimention_one1 - dimention_one0);
    rec.v = (z - dimention_two0) / (dimention_two1 - dimention_two0);
    rec.t = t;
    rec.mat_ptr = matptr;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(1, 0, 0);
    if (dot(r.direction(), rec.normal) > 0) rec.normal = -rec.normal;
    return true;
  }
};

#endif
