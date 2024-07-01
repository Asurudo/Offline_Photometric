#ifndef TRANSFORMATIONH
#define TRANSFORMATIONH

#include "hitable.h"

class translate : public hitable {
 private:
  hitable* hitptr;
  vec3 offset;

 public:
  translate(hitable* hitptr, const vec3& offset)
      : hitptr(hitptr), offset(offset) {}
  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

bool translate::hit(const ray& r, double tmin, double tmax,
                    hit_record& rec) const {
  // 将射线向反方向平移，如果平移后的射线能击中物体，则更改击中点
  // 例如物体本来应该向右平移2个单位，则射线如果向左平移2个单位能击中原来未平移的物体，则默认射线击中了平移后的物体
  ray moved_r(r.origin() - offset, r.direction(), r.time());
  if (hitptr->hit(moved_r, tmin, tmax, rec)) {
    rec.p += offset;
    return true;
  } else
    return false;
}

std::shared_ptr<aabb> translate::getaabb() const {
  std::shared_ptr<aabb> aabbbox = hitptr->getaabb();
  return std::make_shared<aabb>(aabbbox->vec3L + offset,
                                aabbbox->vec3R + offset);
}

class rotate_y : public hitable {
 private:
  double sintheta, costheta;
  hitable* hitptr;
  std::shared_ptr<aabb> aabbboxptr;

 public:
  rotate_y(hitable* hitptr, double angle);
  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

rotate_y::rotate_y(hitable* hitptr, double angle) : hitptr(hitptr) {
  // 转成弧度
  double radians = (PI / 180) * angle;
  sintheta = sin(radians), costheta = cos(radians);
  vec3 vec3L(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 vec3R(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  aabbboxptr = hitptr->getaabb();
  // 计算旋转后的包围盒，遍历八个顶点，计算出八个顶点旋转后的位置
  for (int i = 0; i <= 1; i++)
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++) {
        double orix =
            i * aabbboxptr->vec3R.x() + (1 - i) * aabbboxptr->vec3L.x();
        double oriy =
            j * aabbboxptr->vec3R.y() + (1 - j) * aabbboxptr->vec3L.y();
        double oriz =
            k * aabbboxptr->vec3R.z() + (1 - k) * aabbboxptr->vec3L.z();

        double newx = costheta * orix + sintheta * oriz;
        double newz = -sintheta * orix + costheta * oriz;

        vec3 curp(newx, oriy, newz);
        for (int c = 0; c < 3; c++) {
          vec3R[c] = std::max(vec3R[c], curp[c]);
          vec3L[c] = std::min(vec3L[c], curp[c]);
        }
      }
    
  aabbboxptr = std::make_shared<aabb>(vec3L, vec3R);
}

bool rotate_y::hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const {
  vec3 rayorigin = r.origin();
  vec3 raydirec = r.direction();
  // 将射线旋转
  rayorigin[0] = costheta * r.origin()[0] - sintheta * r.origin()[2];
  rayorigin[2] = sintheta * r.origin()[0] + costheta * r.origin()[2];
  raydirec[0] = costheta * r.direction()[0] - sintheta * r.direction()[2];
  raydirec[2] = sintheta * r.direction()[0] + costheta * r.direction()[2];

  ray rotatedray(rayorigin, raydirec, r.time());
  // 如果旋转后的射线射中，则将击中记录进行旋转
  if(hitptr->hit(rotatedray, tmin, tmax, rec)){
    vec3 point = rec.p;
    vec3 normal = rec.normal;
    point[0] = costheta*rec.p[0] + sintheta*rec.p[2];
    point[2] = -sintheta*rec.p[0] + costheta*rec.p[2];
    normal[0] = costheta * rec.normal[0] + sintheta * rec.normal[2];
    normal[2] = -sintheta * rec.normal[0] + costheta * rec.normal[2];
    rec.p = point;
    rec.normal = normal;
    return true;
  }
  return false;
}

std::shared_ptr<aabb> rotate_y::getaabb() const {
  return std::make_shared<aabb>(aabbboxptr->vec3L, aabbboxptr->vec3R);
}



#endif
