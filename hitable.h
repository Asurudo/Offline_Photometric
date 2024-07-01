#ifndef HITABLEH
#define HITABLEH

#include "ray.h"

// 前向声明
class material;
class aabb;

struct hit_record
{
    // a+tb中的t
    double t;
    // 交点
    vec3 p;
    // 视线与物体相交的点处的法线
    vec3 normal;
    // uv贴图
    double u, v;
    // 材质信息指针
    material* mat_ptr;
};

class hitable
{
    // 纯虚函数
    public:
        virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const = 0;
        // 获得包围盒
        virtual std::shared_ptr<aabb> getaabb() const = 0;
};

#endif