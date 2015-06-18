struct index 
{
    // ...
    typedef int index_type[2];
    const index_type& operator[](int i);
    // ...
};

int main() {
int k = 0;
int i = (index()[k])[1];
}
