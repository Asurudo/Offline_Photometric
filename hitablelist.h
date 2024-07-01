#ifndef HITABLELISTH
#define HITABLELISTH

#include <vector>

#include "bvh.h"
#include "hitable.h"
#include "sphere.h"

class hitable_list {
 public:
  std::vector<std::shared_ptr<hitable>> list;
  std::shared_ptr<hitable> bvhroot;
  hitable_list() {}
  hitable_list(const std::vector<std::shared_ptr<hitable>> l) : list(l) {}
  hitable_list(std::shared_ptr<hitable> bvhroot):bvhroot(bvhroot){}
//   hitable_list& operator=(const hitable_list& world) {
//     list = world.list;
//     bvhroot = world.bvhroot;
//     return *this;
//   }
  virtual bool hitanything(const ray& r, double tmin, double tmax,
                           hit_record& rec) const;
  virtual bool hitanythingbvh(const ray& r, double tmin, double tmax,
                              hit_record& rec) const;
};

bool hitable_list::hitanything(const ray& r, double tmin, double tmax,
                               hit_record& rec) const {
  hit_record tempRec;
  bool hitAnything = false;
  // closestSoFar
  // 保证每次查询的最远距离：如果先遇到近的以后就不会考虑远的，从而不会把远近搞错。
  double closestSoFar = tmax;
  for (const auto& h : list)
    if (h->hit(r, tmin, closestSoFar, tempRec)) {
      hitAnything = true;
      closestSoFar = tempRec.t;
      rec = tempRec;
    }
  return hitAnything;
}

bool hitable_list::hitanythingbvh(const ray& r, double tmin, double tmax,
                                  hit_record& rec) const {
  hit_record tempRec;
  bool hitAnything = false;
  // 直接射向根节点，递归程序会带着往左右子节点遍历
  // 此处tmin和tmax没有意义，因为会计算包围盒的击中时间（在bvh的hit中保证）
  if (bvhroot->hit(r, tmin, tmax, tempRec)) {
    hitAnything = true;
    rec = tempRec;
  }
  return hitAnything;
}

#endif