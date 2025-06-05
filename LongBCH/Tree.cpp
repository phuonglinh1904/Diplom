//
// Created by ADMIN23 on 25.02.2025.
//
#include <memory>
#include <vector>      // Bổ sung cho vector
#include <iostream>    // Bổ sung cho cout
#include "Polynomial.h"
#include "Tree.h"

using namespace std;

Node::Node(const Polynomial& poly) : polynomial(poly), left(nullptr), right(nullptr) {}
shared_ptr<Node> buildSubproductTree(const vector<Polynomial>& polynomials, long start, long end) {
    if (start > end) {
        return nullptr;
    }
    if (start == end) {
        return make_shared<Node>(polynomials[start]);
    }
    long mid = start + (end - start) / 2;
    auto leftChild = buildSubproductTree(polynomials, start, mid);
    auto rightChild = buildSubproductTree(polynomials, mid + 1, end);
    Polynomial product = Polynomial::mul(leftChild->polynomial, rightChild->polynomial);
    auto node = make_shared<Node>(product);
    node->left = leftChild;
    node->right = rightChild;
    return node;
}

void computeRemainders(const Polynomial& h, const shared_ptr<Node>& node, vector<Polynomial>& remainders, const mpz_class &q) {
    if (node->left == nullptr && node->right == nullptr) {
        // Nếu là nút lá, tính phần dư và lưu vào danh sách
        Polynomial remainder = Polynomial::mod(h, node->polynomial, q);
        remainders.push_back(remainder);
    } else {
        Polynomial remainder = Polynomial::mod(h, node->polynomial, q);
        computeRemainders(remainder, node->left, remainders, q);
        computeRemainders(remainder, node->right, remainders, q);
    }
}



void printTree(const shared_ptr<Node>& node, int level) {
    if (!node) return;

    for (int i = 0; i < level; ++i) cout << "    ";
    cout << node->polynomial << endl;
    printTree(node->left, level + 1);
    printTree(node->right, level + 1);
}

