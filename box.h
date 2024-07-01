#ifndef BOXH
#define BOXH

#include <memory>
#include "rectangle.h"
#include "vec3.h"
#include "hitable.h"

class box : public hitable {
 private:
  std::vector<hitable*> reclist;
  // 左下角和右上角坐标
  vec3 pmin, pmax;
  material* matptr;

 public:
  box() {}
  box(const vec3& pmin, const vec3& pmax, material* matptr)
      : pmin(pmin), pmax(pmax), matptr(matptr) {
    reclist.emplace_back(new rectangle_xy(pmin.x(), pmax.x(), pmin.y(),
                                          pmax.y(), pmax.z(), matptr));
    reclist.emplace_back(new rectangle_xy(pmin.x(), pmax.x(), pmin.y(),
                                          pmax.y(), pmin.z(), matptr));
    reclist.emplace_back(new rectangle_xz(pmin.x(), pmax.x(), pmin.z(),
                                          pmax.z(), pmax.y(), matptr));
    reclist.emplace_back(new rectangle_xz(pmin.x(), pmax.x(), pmin.z(),
                                          pmax.z(), pmin.y(), matptr));
    reclist.emplace_back(new rectangle_yz(pmin.y(), pmax.y(), pmin.z(),
                                          pmax.z(), pmax.x(), matptr));
    reclist.emplace_back(new rectangle_yz(pmin.y(), pmax.y(), pmin.z(),
                                          pmax.z(), pmin.x(), matptr));
  }
  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

std::shared_ptr<aabb> box::getaabb() const {
  return std::make_shared<aabb>(pmin, pmax);
}

bool box::hit(const ray& r, double tmin, double tmax, hit_record& rec)const  {
  hit_record tempRec;
  bool hitAnything = false;
  // closestSoFar
  // 保证每次查询的最远距离：如果先遇到近的以后就不会考虑远的，从而不会把远近搞错。
  double closestSoFar = tmax;
  for (const auto& h : reclist)
    if (h->hit(r, tmin, closestSoFar, tempRec)) {
      hitAnything = true;
      closestSoFar = tempRec.t;
      rec = tempRec;
    }
  return hitAnything;
}

#endif
