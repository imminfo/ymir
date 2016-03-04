import os

if __name__ == "__main__":
    for f in os.listdir("./build/test/"):
        if f.find("Test_") != -1:
            print("==================================")
            print("Run test cases from ", f)
            print("==================================")
            print()
            os.system("./build/test/" + f + " " + os.getcwd() + '/test/data/')
            print()
