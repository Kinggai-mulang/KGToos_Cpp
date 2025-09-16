#include<iostream>
#include"FileRead.h"

// const std::string folderPath = "E:\\CDTH\\Data\\TH1217\\Trees";

int main()
{
    //std::string folderpath = "F:\\CDTH\\202406@CDTH@CropperGNSS\\ROVE_20250425_08_THAM03_01";
    //setTHOptfile(folderpath);
    //mergeBASEObsFiles("F:\\CDTH\\202406@CDTH@CropperGNSS\\ROVE_20250425_08_THAM03_01");
    //double XYZ[3];
    //double BLH[3] = { 30.527790303640*D2R,  114.3569040975299*D2R,  81.1906659 };
    //LLH2XYZ(BLH, XYZ);
    if (1) {
        test_incremental_vs_direct();
            
        //test_multi_removal();
    }


    if (0) {
        DecodeNEMA(1);
    }


    if (0) {
        ComparePOS();
    }



    if (0)
    {
        DifPos();
    }

    if(0)
    {
        getRefPos(1);
    }
    
    //
    if (0)
    {
        std::string folderPath;
        std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
        std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
        std::string outputBASEPath = FILEOpt::getParentDir(folderPath) + "\\Base.rtcm3";
        std::string outputROVEPath = FILEOpt::getParentDir(folderPath) + "\\Rove.rtcm3";
        extractRoveFromRAW(folderPath, outputBASEPath, "$RTK_BASE,");
        extractRoveFromRAW(folderPath, outputROVEPath, "$RTK_RAW,");
    }

    if (0) {
        std::string folderPath;
        std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
        std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
        //std::string folderPath = "F:\\CDTH\\2025@Couple\\20250310@Desaysv@rtk@dr@Tunnel\\pos\\location.log";
        std::string outputROVEPath = FILEOpt::getParentDir(folderPath) + "\\location.pos";
        ExtractGNGGA(folderPath, outputROVEPath);
    }


    if (0) {

        std::string folderPath;
        std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
        std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
        extraDaoyYuan(folderPath);
    }
    //FILEOpt::changeFileExtension(outputPath, ".rtcm3");

    ////std::cout << std::endl << "成功解码" << std::endl;
    //}
    

    //uint8_t a = 0;
    //std::string folderPath;
    //std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
    //std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
    //std::string outputPath = FILEOpt::getParentDir(folderPath)+"\\ROVE.rtcm3";
    ////FILEOpt::changeFileExtension(outputPath, ".rtcm3");
    //extractGNGGA(folderPath, outputPath);
    //std::cout << "是否为单边场景？   如果为单边设置为输入1，否则输入0" << std::endl;
    //std::cin >> a;
    //
    //auto Pfile = FILEOpt::traversePdir(folderPath);
    //for (auto sfile = Pfile.begin(); sfile != Pfile.end(); sfile++)
    //{
    //    if(a!='1')FILEOpt::changeFileExtension(*sfile, ".rtcm3");
    //    if(a=='1')renameTH_OneWallsDATA(*sfile);
    //}

    
    //RenameObsEphFILE_ROVE(folderPath);
    //PDirTHobs(folderPath);
    //folderPath = FILEOpt::getParentDir(folderPath, 1);
    //setTHOptfile(folderPath);

    if (0) {
        std::string folderPath;
        std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
        std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
        setTHOptfile(folderPath);
    }

    if (0) {
        std::string folderPath;
        std::cout << "请输入文件夹路径:    (原始观测数据的文件夹路径，如THdebug)" << std::endl;
        std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
        //renameTHObsFile(folderPath);
        //RenameObsEphFILE_ROVE(folderPath);
        //PDirTHobs(folderPath);
        //folderPath = FILEOpt::getParentDir(folderPath, 1);
        setTHOptfile(folderPath);
    }


        

    //std::map<std::string, std::vector<std::string>> file;
    //std::string output;

    //file = FILEOpt::traverseGPdir(folderPath);
    //auto Pfile = FILEOpt::traversePdir(folderPath);
    //////batchExecute(Pfile);
    //for (auto sfile =Pfile.begin(); sfile !=Pfile.end(); sfile++)
    //{
    //    //renameTHObsFile(sfile[0]);
    //    //FILEOpt::changeFileExtension(*sfile, ".rtcm3");
    //    renameTH_OneWallsDATA(*sfile);
    //}

    //for (const auto& pair : file) {
    //    // pair.first 是键 (std::string)
    //    // pair.second 是值 (std::vector<std::string>)

    //    std::cout << "Folder: " << pair.first << std::endl;
    //    for (const auto& file : pair.second) {
    //        std::cout << "  File: " << file << std::endl;
    //        //std::string outFolder = FILEOpt::getParentDir(file,2);
    //        //outFolder = FILEOpt::CreatFolder(outFolder, "V313Degug");
    //        //FILEOpt::copyFileToFolder(file, outFolder);
    //        //output = FILEOpt::changeStrExtension(file,".pos");
    //        //output = FILEOpt::getParentFileName(output, 0);  //这里仅仅得到文件名
    //        //output = outFolder +"\\"+ output;
    //        ////copyBinFile(file, output);
    //        //extractGNGGA(file, output);
    //        //FILEOpt::changeFileExtension(file, ".rtcm3");
    //        renameTHObsFile(file);
    //    }
    //}
    //std::cout << "遍历文件夹数量   " << FILEOpt::Foldercount << std::endl;
    //std::cout << "遍历文件数量    " << FILEOpt::Filecount << std::endl;


    //system("pause");

    return 0;
}





