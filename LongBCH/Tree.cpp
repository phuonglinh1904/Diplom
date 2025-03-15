//
// Created by ADMIN23 on 25.02.2025.
//
#include <memory>
#include <vector>      // Bổ sung cho vector
#include <iostream>    // Bổ sung cho cout
#include "Polynomial.h"
#include "Tree.h"

using namespace std;

// Định nghĩa constructor cho Node
Node::Node(const Polynomial& poly) : polynomial(poly), left(nullptr), right(nullptr) {}

// Hàm xây dựng cây tích con
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

    // Tính tích của hai đa thức con
    Polynomial product = Polynomial::mul(leftChild->polynomial, rightChild->polynomial);
    auto node = make_shared<Node>(product);
    node->left = leftChild;
    node->right = rightChild;
    return node;
}

// Hàm tính phần dư của h_i(x) với g trên cây subproduct
void computeRemainders(const Polynomial& h, const shared_ptr<Node>& node, vector<Polynomial>& remainders, const mpz_class &q) {
    if (node->left == nullptr && node->right == nullptr) {
        // Nếu là nút lá, tính phần dư và lưu vào danh sách
        Polynomial remainder = Polynomial::mod(h, node->polynomial, q);
        remainders.push_back(remainder);
    } else {
        // Nếu là nút trong, tính phần dư và đệ quy xuống các nút con
        Polynomial remainder = Polynomial::mod(h, node->polynomial, q);
        computeRemainders(remainder, node->left, remainders, q);
        computeRemainders(remainder, node->right, remainders, q);
    }
}



void printTree(const shared_ptr<Node>& node, int level) {
    if (!node) return;

    // In nút hiện tại
    for (int i = 0; i < level; ++i) cout << "    ";
    cout << node->polynomial << endl;

    // Đệ quy in cây con trái và phải
    printTree(node->left, level + 1);
    printTree(node->right, level + 1);
}

