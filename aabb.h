#ifndef AABBH
#define AABBH

#include "hitable.h"

// 包围盒
class aabb : public hitable {
  // 包围盒的左右上下前后边界
 public:
  vec3 vec3L, vec3R;
  aabb(){
    vec3L = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
    vec3R = vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  }

  // 将两个小包围盒组建成更大的包围盒，左边界取两个边界的最小值，右边界取两个边界的最大值
  aabb(const aabb& box1, const aabb& box2) {
    vec3L = vec3(std::min(box1.vec3L[0], box2.vec3L[0]),
                 std::min(box1.vec3L[1], box2.vec3L[1]),
                 std::min(box1.vec3L[2], box2.vec3L[2]));
    vec3R = vec3(std::max(box1.vec3R[0], box2.vec3R[0]),
                 std::max(box1.vec3R[1], box2.vec3R[1]),
                 std::max(box1.vec3R[2], box2.vec3R[2]));
  }
  aabb(const vec3& l, const vec3& r) : vec3L(l), vec3R(r) {}

  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

// 判断r有没有击中这个包围盒
// 线段相交问题，xyz每个维度都对应t0到t1的参数
// 和盒子相交，要求存在一个时间（参数）t，属于三个维度的[t0,t1]
bool aabb::hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
  for (int i = 0; i < 3; i++) {
    double invD = 1.0 / r.direction()[i];
    double t0 = (vec3L[i] - r.origin()[i]) * invD;
    double t1 = (vec3R[i] - r.origin()[i]) * invD;
    if (invD < 0.0f) std::swap(t0, t1);
    // 统计最大的左端点与最小的右端点
    tmin = std::max(t0, tmin);
    tmax = std::min(t1, tmax);
    // 如果最小的右端点比最大的左端点小，说明至少有两条线段没有重叠部分
    if (tmax <= tmin) return false;
  }
  // 进入包围盒的时间
  rec.t = tmin;
  return true;
}

std::shared_ptr<aabb> aabb::getaabb() const {
  return std::make_shared<aabb>(vec3L, vec3R);
}

#endif