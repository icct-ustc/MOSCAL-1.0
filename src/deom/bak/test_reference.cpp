#include <iostream>
#include <vector>

using namespace std;

int main() {
  vector<int> a(10);
  vector<int> b(10);
  vector<int> &c = a;
  c[0] = 1;
  c = b;
  c[0] = 10;
  for (auto i : a) {
    cout << i << endl;
  }
  for (auto i : b) {
    cout << i << endl;
  }
  return 0;
}
