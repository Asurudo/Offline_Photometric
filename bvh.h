#ifndef BVHH
#define BVHH

#include <cassert>
#include <memory>

#include "aabb.h"
#include "hitable.h"
#include "hitablelist.h"
#include "jyorand.h"

// 前向声明
extern Rand jyorandengine;

class bvh_node : public hitable {
 private:
  // BVH树的左右子树
  std::shared_ptr<hitable> leftptr;
  std::shared_ptr<hitable> rightptr;
  // 该BVH节点对应的包围盒
  std::shared_ptr<aabb> box;

  // 根据X,Y,Z轴的左边界大小进行排序的，三个不同的比较函数
  static bool box_compare(const std::shared_ptr<hitable> a,
                          const std::shared_ptr<hitable> b, int axis_index) {
    return a->getaabb()->vec3L[axis_index] < b->getaabb()->vec3L[axis_index];
  }

  static bool box_x_compare(const std::shared_ptr<hitable> a,
                            const std::shared_ptr<hitable> b) {
    return box_compare(a, b, 0);
  }

  static bool box_y_compare(const std::shared_ptr<hitable> a,
                            const std::shared_ptr<hitable> b) {
    return box_compare(a, b, 1);
  }

  static bool box_z_compare(const std::shared_ptr<hitable> a,
                            const std::shared_ptr<hitable> b) {
    return box_compare(a, b, 2);
  }

 public:
  bvh_node() {}
  // 传入一个可碰撞物体的列表，返回root，即BVH树的最大包围盒节点
  bvh_node(const std::vector<std::shared_ptr<hitable>>& list,
           std::shared_ptr<hitable>& root) {
    if (list.size() == 0) assert(0 == 1);
    root = getroot(list, 0, list.size() - 1);
  }
  
  // 建立BVH树的算法
  std::shared_ptr<hitable> getroot(
      const std::vector<std::shared_ptr<hitable>>& list, size_t start,
      size_t end) {
    // 将列表复制一份以便排序
    std::vector<std::shared_ptr<hitable>> objlist(list);

    // 随机选择一个轴，将其作为基准对物体进行排序
    int axis = jyorandengine.jyoRandGetInteger(0, 2);
    auto comparator = (axis == 0)
                          ? box_x_compare
                          : ((axis == 1) ? box_y_compare : box_z_compare);
    size_t objectsize = end - start + 1;

    // 目前的节点
    std::shared_ptr<bvh_node> rootptr = std::make_shared<bvh_node>();
    if (objectsize == 1) {
      // 如果只有一个物体，则左右指针都是sphere类型
      rootptr->leftptr = rootptr->rightptr = objlist[start];
    } else {
      // 否则进行排序后，递归下去
      std::sort(objlist.begin() + start, objlist.begin() + end + 1, comparator);
      size_t mid = (start + end) / 2;
      rootptr->leftptr = getroot(objlist, start, mid);
      rootptr->rightptr = getroot(objlist, mid + 1, end);
    }

    // 将左右的小包围盒组成更大的包围盒，赋值给目前的节点
    rootptr->box = std::make_shared<aabb>(*(rootptr->leftptr->getaabb()),
                                          *(rootptr->rightptr->getaabb()));
    return rootptr;
  }

  virtual bool hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const override;
  virtual std::shared_ptr<aabb> getaabb() const override;
};

// 判断是否穿过该节点对应包围盒
bool bvh_node::hit(const ray& r, double tmin, double tmax,
                   hit_record& rec) const {
  hit_record rc;
  
  // 取出包围盒进行aabb的hit判断，如果没穿过则直接返回false
  if (!getaabb()->hit(r, tmin, tmax, rc)) return false;
  
  // 判断是否穿过左节点对应包围盒
  bool hit_left = leftptr->hit(r, tmin, tmax, rec);
  // 判断是否穿过右节点对应包围盒，注意，若穿过了左节点，则进入左节点包围盒的时间t则为右节点的tmax
  bool hit_right = rightptr->hit(r, tmin, hit_left ? rec.t : tmax, rec);
  
  // 至少需要经过一个子节点
  return hit_left || hit_right;
}

// bvh节点本身不判断碰撞，返回包围盒
std::shared_ptr<aabb> bvh_node::getaabb() const {
  return std::make_shared<aabb>(box->vec3L, box->vec3R);
}

#endif