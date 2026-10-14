#pragma once
#include "library/geometry/pointbase.hpp"
// 偏角ソートに対する argsort
// -pi~piの順で同じ角度の場合は順番は決まっていない
template <typename T>
vector<int> angle_sort(vector<PointBase<T>>& P) {
    vector<int> lower, origin, upper;
    const PointBase<T> O = {0, 0};
    for(int i = 0;i < P.size();i++) {
        if (P[i] == O) origin.emplace_back(i);
        else if ((P[i].y < 0) || (P[i].y == 0 && P[i].x > 0)) lower.emplace_back(i);
        else upper.emplace_back(i);
    }
    sort(all(lower), [&](auto& i, auto& j) { return cross(P[i],P[j]) > 0; });
    sort(all(upper), [&](auto& i, auto& j) { return cross(P[i],P[j]) > 0; });
    auto& I = lower;
    I.insert(I.end(), origin.begin(),origin.end());
    I.insert(I.end(), upper.begin(),upper.end());
    return I;
}


template <typename T>
vector<int> angle_sort(vector<pair<T, T>>& P) {
    vector<PointBase<T>> tmp(P.size());
    for(int i = 0;i < P.size();i++) tmp[i] = PointBase<T>(P[i]);
    return angle_sort<T>(tmp);
}