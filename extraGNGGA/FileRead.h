#pragma once
#include <filesystem>
#include <cstring>
#include<fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <functional>
#include"BaseGT.h"
#include"CMatrix.h"
constexpr int MAXSIZE = 1023;   //最大字符长度
//constexpr double PI = 3.1415926535897932384626433832795;		///< pi

//针对int类型返回的数据进行说明
enum class ErrorCode {
    SUCCESS = 0,
    FILE_NOT_FOUND = 1,
    FILE_SIZE_ZERO = 2,
    FOLDER_NOT_FOUND=3,
    FOLDER_SIZE_ZERO=4,
    INVALID_DATA = -1,
    UNSUPPORTED_TYPE = -999
};


class cStrOpt
{
public:
    static double tod(const char* str, const int nPos=0, const int nCount=MAXSIZE);
    static int toi(const char* str, const int nPos=0, const int nCount=MAXSIZE);
    static size_t toui(const char *str, const int nPos = 0, const int nCount = MAXSIZE);

    static void RemoveFileExtension(char* path);
    static void ReplaceFilePath(char* FilePath);
    static void xstrmid(const char* src, const int nPos, const int nCount, char* dst);
};



//主要用来实现文件的各种常用操作!
class FILEOpt
{
public:
    static std::string  lineIdentifier;             //行内标识符号
    static std::string  FolderIdntifier;
    static std::string  FileIdentifier /*= "."*/;    // 标识符,文件中得有LC29H
    static std::string  Filesuffix /*= ".pos"*/;     // 后缀表示符，处理后缀名.pos文件
   
public:
    static  int         Filecount;               // 文件计数，记录处理了多少个文件
    static  int         Foldercount;             // 文件夹计数，记录处理了多少个文件夹

    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    static std::vector<std::string> traversePdir(const std::string& folderPath);
    static std::map<std::string, std::vector<std::string>> traverseGPdir(const std::string& rootFolderPath);  // 遍历更高级别目录
    static bool traverseALLFolder(const std::string& folderPath, const std::function<void(const std::string&)>& fileOperation);   //遍历所有文件目录

    static std::string getParentDir(const std::string& filePath, const short n = 1);
    static std::string getParentFileName(const std::string& filePath, const short n = 1);
    static bool trim(std::string& value, bool frist = true, bool last = true);

    // 类相关
    static int getFileCount();
    static int getFolderCount();

    static bool FileExists(const std::string& filePath);  //  判断文件是否存在
    static bool SkipLines(std::ifstream& file, const int linesToSkip);      // 文件跳过多少行

    static std::string CreatFolder(const std::string& path, const std::string Foldername);   // 生成文件夹
    static void changeFileExtension(const std::string& filePath, const std::string& newExtension);  //改变文件后缀名
    static bool copyFileToFolder(const std::string &fpath,const  std::string &oFolder);

    static uintmax_t getFileSize(const std::string& file);

    static std::string changeStrExtension(const std::string& filePath, const std::string& extension);   //改变string类型的后缀名.pos
    static void copyBinFile(const std::string& inputFilePath, const std::string& outputFilePath);   // 将二进制文件复制到一个新文件
    
    //
};

bool DecodeNEMA(int PosType=0);


void LLH2XYZ(const double LLH[3], double XYZ[3]);
double DegMinSec2D(const double& X);
void XYZ2LLH(const double XYZ[3], double LLH[3]);
double meanWithoutOutliers(const std::vector<double>& data);



void getRefPos(int PosType=0);




void extraDaoyYuan(const std::string& filePath);
void extractRoveFromRAW(const std::string& inputFilePath, const std::string& outputFilePath, const char SimbolIdt[]);
bool ExtractGNGGA(const std::string& inputFile, const std::string& outputFile);
bool renameTHObsFile(const std::string& oldFilePath);


bool PDirTHobs(std::string folderPath);



void renameTH_OneWallsDATA(std::string& filepath);

void setTHOptfile(const std::string& folderPath);

void mergeBASEObsFiles(const std::filesystem::path& dir);

void RenameObsEphFILE_ROVE(const std::string& folderPath);

void DifPos(/*std::string &InputIPS,std::string &InputCouple*/);


void ComparePOS();



// Test
void test_incremental_vs_direct();
void generate_positive_definite_matrix(int n, double* matrix);
void test_multi_removal();