#include <boost/multiprecision/cpp_int.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
namespace mp = boost::multiprecision;

// Elliptic curve parameters (secp256k1)
const mp::cpp_int P("115792089237316195423570985008687907853269984665640564039457584007908834671663");
const mp::cpp_int N("115792089237316195423570985008687907852837564279074904382605163141518161494337");
const mp::cpp_int A = 0;
const mp::cpp_int B = 7;
const mp::cpp_int Gx("55066263022277343669578718895168534326250603453777594175500187360389116729240");
const mp::cpp_int Gy("32670510020758816978083085130507043184471273380659243275938904335757337482424");

using Point = std::tuple<mp::cpp_int, mp::cpp_int, mp::cpp_int>;

mp::cpp_int inv(mp::cpp_int a, mp::cpp_int n) {
    if (a == 0) return 0;
    mp::cpp_int lm = 1, hm = 0;
    mp::cpp_int low = a % n, high = n;
    while (low > 1) {
        mp::cpp_int r = high / low;
        mp::cpp_int nm = hm - lm * r;
        mp::cpp_int new_ = high - low * r;
        hm = lm;
        high = low;
        lm = nm;
        low = new_;
    }
    return lm % n;
}

Point to_jacobian(Point p) {
    return make_tuple(get<0>(p), get<1>(p), 1);
}

Point jacobian_double(Point p) {
    if (get<1>(p) == 0) return make_tuple(0, 0, 1);
    mp::cpp_int ysq = (get<1>(p) * get<1>(p)) % P;
    mp::cpp_int S = (4 * get<0>(p) * ysq) % P;
    mp::cpp_int M = (3 * get<0>(p) * get<0>(p) + A) % P;
    mp::cpp_int nx = (M * M - 2 * S) % P;
    mp::cpp_int ny = (M * (S - nx) - 8 * ysq * ysq) % P;
    return make_tuple(nx, ny, 2 * get<1>(p) * get<2>(p) % P);
}

Point jacobian_add(Point p, Point q) {
    if (get<1>(p) == 0) return q;
    if (get<1>(q) == 0) return p;
    mp::cpp_int U1 = (get<0>(p) * get<2>(q) * get<2>(q)) % P;
    mp::cpp_int U2 = (get<0>(q) * get<2>(p) * get<2>(p)) % P;
    mp::cpp_int S1 = (get<1>(p) * get<2>(q) * get<2>(q) * get<2>(q)) % P;
    mp::cpp_int S2 = (get<1>(q) * get<2>(p) * get<2>(p) * get<2>(p)) % P;
    if (U1 == U2) {
        if (S1 != S2) return make_tuple(0, 0, 1);
        return jacobian_double(p);
    }
    mp::cpp_int H = U2 - U1;
    mp::cpp_int R = S2 - S1;
    mp::cpp_int H2 = (H * H) % P;
    mp::cpp_int H3 = (H * H2) % P;
    mp::cpp_int U1H2 = (U1 * H2) % P;
    mp::cpp_int nx = (R * R - H3 - 2 * U1H2) % P;
    mp::cpp_int ny = (R * (U1H2 - nx) - S1 * H3) % P;
    mp::cpp_int nz = (H * get<2>(p) * get<2>(q)) % P;
    return make_tuple(nx, ny, nz);
}

Point from_jacobian(Point p) {
    mp::cpp_int z = inv(get<2>(p), P);
    return make_tuple((get<0>(p) * z * z) % P, (get<1>(p) * z * z * z) % P, 1);
}

Point jacobian_multiply(Point a, mp::cpp_int n) {
    if (get<1>(a) == 0 || n == 0) return make_tuple(0, 0, 1);
    if (n == 1) return a;
    if (n < 0 || n >= N) return jacobian_multiply(a, n % N);
    if (n % 2 == 0) return jacobian_double(jacobian_multiply(a, n / 2));
    return jacobian_add(jacobian_double(jacobian_multiply(a, n / 2)), a);
}

Point multiply(Point a, mp::cpp_int n) {
    return from_jacobian(jacobian_multiply(to_jacobian(a), n));
}
std::string cppIntToHexString(const mp::cpp_int& value) {
    std::stringstream ss;
    ss << std::hex << value;
    return ss.str();
}

int main() {
  
  std::string hexString = "c398a7e540bc209884c113a0b8b018c1a112367f0d5ce0023c162cc8e5f8ebbe";

  // Convert hexadecimal string to decimal
  mp::cpp_int decimalValue;
  std::stringstream hexStream;
  hexStream << std::hex << hexString;
  hexStream >> decimalValue;

  // Compute public key
  auto Pub = multiply(make_tuple(Gx, Gy, 1), decimalValue);

  mp::cpp_int pub_x = get<0>(Pub);
  mp::cpp_int pub_y = get<1>(Pub);

  std::string Pub_F = cppIntToHexString(pub_x) + cppIntToHexString(pub_y);
  std::cout << Pub_F << std::endl,

  return 0;
}
