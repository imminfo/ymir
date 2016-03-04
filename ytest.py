import os

if __name__ == "__main__":
    for f in os.listdir("./build/test/"):
        if f.find("Test_") != -1:
            os.system(f + os.getcwd() + '/test/data/')
