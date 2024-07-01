#ifndef SMOKEH
#define SMOKEH

#include "hitable.h"
#include "material.h"

extern Rand jyorandengine;

// 将一个物体烟雾化
class smoke : public hitable {
 private:
  hitable* hitptr;
  double density;
  isotropic* isotropicptr;

 public:
  smoke(hitable* hitptr, double density, texture* textureptr)
      : hitptr(hitptr), density(density) {
    isotropicptr = new isotropic(textureptr);
  }
  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

std::shared_ptr<aabb> smoke::getaabb() const { return hitptr->getaabb(); }

bool smoke::hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
  hit_record rec1, rec2;
  // 射两次，第一次找入口，第二次找出口
  // 为什么是从-FLT_MAX开始而不是从0开始，因为在烟雾内，的确可能射回去（反方向）
  if (hitptr->hit(r, -FLT_MAX, FLT_MAX, rec1) &&
      hitptr->hit(r, rec1.t + 0.0001, FLT_MAX, rec2)) {
    
    // 确保入口的t小于出口的t，且在tmin和tmax范围内
    rec1.t = std::max(tmin, rec1.t);
    rec2.t = std::min(tmax, rec2.t);
    if (rec1.t >= rec2.t) return false;
    rec1.t = std::max(rec1.t, 0.0);

    // 光线穿过烟雾（从入口到出口）的所需要的总距离
    double distance_inside_boundary =
        (rec2.t - rec1.t) * r.direction().length();
    
    // 光线在烟雾内部某点的碰撞折射所需要的距离
    // 密度越高，碰撞概率越大
    double hit_distance =
        -(1 / density) * log(jyorandengine.jyoRandGetReal<double>(0, 1));

    // 碰撞折射距离小于总距离，表明发生了碰撞
    if (hit_distance < distance_inside_boundary) {
      rec.t = rec1.t + hit_distance / r.direction().length();
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = vec3(jyorandengine.jyoRandGetReal<double>(0, 1),
                        jyorandengine.jyoRandGetReal<double>(0, 1),
                        jyorandengine.jyoRandGetReal<double>(0, 1));
      rec.mat_ptr = isotropicptr;
      return true;
    }
  }
  return false;
}

#endif
