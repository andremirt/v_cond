struct proxy {
    operator int(); // int function
    operator double(); // double function
    // proxy(arguments);
    // arguments &arguments_;
};

proxy function() {
    return proxy();
}

int main() {
    int v = function();
    double u  = function();
}
