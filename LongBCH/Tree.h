//
// Created by ADMIN23 on 25.02.2025.
//

#ifndef LONGBCH_TREE_H
#define LONGBCH_TREE_H
#include <memory>
#include <vector>
#include <iostream>
#include "Polynomial.h"
using namespace std;
// Định nghĩa cấu trúc Node của cây subproduct
struct Node {
    Polynomial polynomial; // Đa thức trong nút
    shared_ptr<Node> left; // Con trái
    shared_ptr<Node> right; // Con phải

    // Constructor cho Node
    Node(const Polynomial& poly);
};

// Hàm xây dựng cây tích con
shared_ptr<Node> buildSubproductTree(const vector<Polynomial>& polynomials, long start, long end);
void computeRemainders(const Polynomial& h, const shared_ptr<Node>& node, vector<Polynomial>& remainders, const mpz_class &q);

// Hàm in cây subproduct theo thứ tự Inorder
void printTree(const shared_ptr<Node>& node, int level = 0);
#endif //LONGBCH_TREE_H
