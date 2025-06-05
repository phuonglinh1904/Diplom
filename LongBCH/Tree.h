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
struct Node {
    Polynomial polynomial;
    shared_ptr<Node> left;
    shared_ptr<Node> right;

    Node(const Polynomial& poly);
};

shared_ptr<Node> buildSubproductTree(const vector<Polynomial>& polynomials, long start, long end);
void computeRemainders(const Polynomial& h, const shared_ptr<Node>& node, vector<Polynomial>& remainders, const mpz_class &q);

void printTree(const shared_ptr<Node>& node, int level = 0);
#endif //LONGBCH_TREE_H
