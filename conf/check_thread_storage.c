#if defined(_MSC_VER)
#define __thread __declspec(thread)
#endif

__thread int x;

int main(int argc, char **argv) {
  return x;
}
