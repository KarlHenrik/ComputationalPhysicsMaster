#include <iostream>
#include <cmath>

using namespace std;

class Book {
    private:
        string title;
        string author;
        int pages;
    public:
        Book(string atitle, string aauthor, int apages) {
            title = atitle;
            author = aauthor;
            pages = apages;
        }
        int getPages() {
            return pages;
        }
};

class Novel : public Book {

    public:
        using Book::Book; //this is nicer, but not sure how any of this works
        //Novel(string atitle, string aauthor, int apages):Book(atitle, aauthor, apages) {
        //
        //}
        int getPages() {
            return 100;
        }
};

int main(int argc, char const *argv[]) {

    Novel book1("Harry", "Me", 500);

    cout << book1.getPages();

    return 0;
}
