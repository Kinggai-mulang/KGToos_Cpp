#include"FileRead.h"
#include <iostream>
#include <thread>
#include <regex>






// 静态成员变量，需要在cpp文件夹下定义
int FILEOpt::Foldercount = 0;
int FILEOpt::Filecount = 0;
std::string FILEOpt::Filesuffix = ".";				//文件后缀设置
std::string FILEOpt::FolderIdntifier = "THdebug";			//寻找只有@符号的文件夹
std::string FILEOpt::FileIdentifier = "";			//文件标识符号
std::string FILEOpt::lineIdentifier = "$GNGGA";			//文件行标识符号



static int skipRinexHead(std::ifstream& file)
{
	std::string line;
	int i = 0;
	if (file.is_open()) {
		while (1)
		{
			i++;
			std::getline(file, line);
			if (line.find("END OF HEADER"))return true;
			if (i > 1000)return false;
		}
	}
}

double cStrOpt::tod(const char* str, const int nPos, const int nCount)
{
	if (!str || nPos < 0 || nCount <= 0) {
		return 0.0; // 转化失败
	}

	// 动态计算 nCount，如果为 MAXSIZE，则使用字符串的实际长度
	int length = 0;
	while (str[length] != '\0') {
		++length;
	}
	int count = (nCount == MAXSIZE) ? length - nPos : nCount;
	if (count <= 0 || nPos >= length) {
		return 0.0;
	}

	double result = 0.0;
	double fraction = 0.1; // 用于小数部分的解析
	bool isNegative = false;
	bool isFraction = false;
	int exponent = 0;      // 指数部分
	bool expNegative = false;
	bool hasExponent = false;

	// 遍历从指定位置开始的子字符串
	for (int i = nPos; i < nPos + count && str[i] != '\0'; ++i) {
		char c = str[i];
		if (c == '-' && i == nPos) { // 负号必须在子字符串的第一个字符
			isNegative = true;
		}
		else if (c == '.' && !hasExponent) { // 小数点
			isFraction = true;
		}
		else if ((c == 'e' || c == 'E') && !hasExponent) { // 科学计数法标志
			hasExponent = true;
			if (str[i + 1] == '-') {
				expNegative = true;
				++i; // 跳过负号
			}
			else if (str[i + 1] == '+') {
				++i; // 跳过正号
			}
		}
		else if (c >= '0' && c <= '9') {
			if (!hasExponent) {
				if (!isFraction) {
					result = result * 10 + (c - '0');
				}
				else {
					result += (c - '0') * fraction;
					fraction *= 0.1;
				}
			}
			else {
				exponent = exponent * 10 + (c - '0');
			}
		}
		else {
			// 非法字符，停止解析
			continue;
		}
	}

	// 应用指数部分
	if (hasExponent) {
		double expValue = pow(10.0, expNegative ? -exponent : exponent);
		result *= expValue;
	}

	return isNegative ? -result : result;
}


int cStrOpt::toi(const char* str, const int nPos, const int nCount)
{
	return (int)tod(str, nPos, nCount);
}

size_t cStrOpt::toui(const char* str, const int nPos, const int nCount)
{
	return fabs(cStrOpt::toi(str,nPos,nCount));
}



void cStrOpt::RemoveFileExtension(char* path)
{
	char* dot = strrchr(path, '.');  // 查找最后一个点的位置
	if (dot != NULL) {
		*dot = '\0';  // 将最后一个点替换为字符串终止符
	}
}

void cStrOpt::ReplaceFilePath(char* FilePath)
{
	if (!FilePath) return;

	int n = (int)strlen(FilePath);
	for (int i = 0; i < n; i++)
	{
		if (FilePath[i] == '/') FilePath[i] = '\\';
	}

	// Remove trailing spaces by checking each character directly
	while (n > 0 && (FilePath[n - 1] == ' ' || FilePath[n - 1] == '\t' ||
		FilePath[n - 1] == '\n' || FilePath[n - 1] == '\v' ||
		FilePath[n - 1] == '\f' || FilePath[n - 1] == '\r')) {
		FilePath[n - 1] = '\0';
		n--;
	}
}

/**
* @brief        读取复制字符串
* @param[in]    src       char    源字符串,如果末位有'\n',自动除去
* @param[in]    nPos      int     读取的初始位置(从0开始)
* @param[in]    nCount    int     读取的字符串数
* @param[out]   dst       char    目标字符串, 末尾已用\0填充
* @return       void
* @note         无
* @internals    无
*/
void cStrOpt::xstrmid(const char* src, const int nPos, const int nCount, char* dst)
{
	int		i;
	const char* str;
	char	c;

	str = src + nPos;

	for (i = 0; i < nCount; i++)
	{
		c = *(str + i);
		if (c)
		{
			*(dst + i) = c;
		}
		else
		{
			// 除去末尾的'\n'
			if (dst[i - 1] == '\n') dst[i - 1] = '\0';
			*(dst + i) = '\0';
			break;
		}
	}

	*(dst + nCount) = '\0';
}


/**
 * @file         "FIleRead.cpp"
 * @brief        分割字符串
 * @param[in]    str
 * @param[in]    delimiter
 * @return       std::vector<std::string>
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
std::vector<std::string> FILEOpt::splitString(const std::string& str, char delimiter)
{
	std::vector<std::string> elements;
	std::stringstream ss(str);
	std::string item;
	while (std::getline(ss, item, delimiter)) {
		trim(item);
		if (item != "")elements.push_back(item);
	}
	return elements;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        遍历父级文件夹parent，返回满足条件的所有文件的文件路径
 * @param[in]    folderPath
 * @return       std::vector<std::string>
 * @note         Filesuffix=="."表示不对后缀名进行限制
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
std::vector<std::string> FILEOpt::traversePdir(const std::string& folderPath)
{
	std::vector<std::string>allfile;
	for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			const std::filesystem::path& filePath = entry.path();
			if (filePath.extension()==Filesuffix && filePath.filename().string().find(FileIdentifier) != std::string::npos || Filesuffix==".") {
				allfile.push_back(filePath.string());
				Filecount++;
			}
		}
	}
	return allfile;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        遍历上上级别目录，grandpa目录
 * @param[in]    rootFolderPath
 * @return       std::map<std::string,std::vector<std::string>>
 * @note         
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
std::map<std::string, std::vector<std::string>> FILEOpt::traverseGPdir(const std::string& rootFolderPath)
{
	std::map<std::string, std::vector<std::string>> folderFileMap;
	std::filesystem::path rootF(rootFolderPath);
	std::string asda = rootF.filename().string();
	// 遍历 rootFolderPath 下的每个子文件夹
	for (const auto& entry : std::filesystem::directory_iterator(rootFolderPath)) {
		if (entry.is_directory()) {  // 只处理文件夹
			const std::string folderName = entry.path().filename().string();
			if (folderName.find(FolderIdntifier) == std::string::npos)continue;
			std::vector<std::string> filesInFolder = traversePdir(entry.path().string());  // 获取该文件夹下所有满足条件的文件

			if (!filesInFolder.empty()&& folderName.find(FolderIdntifier)!=std::string::npos) {
				folderFileMap[folderName] = filesInFolder;  // 将文件夹名作为 key，文件路径列表作为 value 存入 map
				Foldercount++;
			}
		}
	}
	return folderFileMap;
}


/**
 * @file         "FIleRead.cpp"
 * @brief        获得某个文件路径的上n级文件路径
 * @param[in]    filePath
 * @param[in]    level
 * @return       std::string
 * @note
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
std::string FILEOpt::getParentDir(const std::string& filePath, const short n)
{
	if (n < 1)return "err";
	// 将路径转换为 std::filesystem::path 类型
	std::filesystem::path path(filePath);
	
	std::filesystem::path parentDir = path.parent_path();
	for (size_t i = 0; i < n-1; i++)
	{
		parentDir = parentDir.parent_path();
	}
	return parentDir.string();
}


/**
 * @file         "FIleRead.cpp"
 * @brief        获得上级目录的名字,n为0的时候，则获取当前文件名
 * @param[in]    filePath
 * @param[in]    n
 * @return       std::string
 * @note         
 * @par History
 *               2024/12/19, Ken Ga  new \n
 */
std::string FILEOpt::getParentFileName(const std::string& filePath, const short n)
{
	if (n < 0)return "err";
	// 将路径转换为 std::filesystem::path 类型
	std::filesystem::path path(filePath);
	if (n == 0)return path.filename().string();
	std::filesystem::path parentDir = path.parent_path();
	for (size_t i = 0; i < n - 1; i++)
	{
		parentDir = parentDir.parent_path();
	}
	return parentDir.filename().string();
}


/**
 * @file         "FIleRead.cpp"
 * @brief        去除前导空白字符，后导空白字符  类似 "   asdad    "->"asdad"
 * @param[in]    value
 * @param[in]    frist
 * @param[in]    last
 * @return       bool
 * @note         如果传入的字符全为空，则会返回""字符，比如传入"       "，则返回"";
 * @par History
 *               2024/12/11, Ken Ga  new \n
 */
bool FILEOpt::trim(std::string& value, bool frist, bool last)
{
	if (frist)   // 去除前导空白字符
	{
		size_t first_non_space = value.find_first_not_of(" \t");
		if (first_non_space != std::string::npos) {
			value = value.substr(first_non_space);  // 去除前导空白字符
		}
		else
			return false;
	}
	if (last)  // 去除后导空白字符
	{
		size_t last_non_space = value.find_last_not_of(" \t");
		if (last_non_space != std::string::npos) {
			value = value.substr(0, last_non_space + 1);  // 去除后导空白字符
		}
		else
			return false;
	}
	return true;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        返回批处理中已经计算的文件数量
 * @return       int
 * @note         
 * @par History
 *               2024/12/13, Ken Ga  new \n
 */
int FILEOpt::getFileCount()
{
	return Filecount;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        返回批处理中已经计算的文件夹数量
 * @return       int
 * @note         
 * @par History
 *               2024/12/13, Ken Ga  new \n
 */
int FILEOpt::getFolderCount()
{
	return Foldercount;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        文件是否存在
 * @param[in]    filePath
 * @return       bool
 * @note         
 * @par History
 *               2024/12/12, Ken Ga  new \n
 */
bool FILEOpt::FileExists(const std::string& filePath)
{
	std::ifstream file(filePath);
	if (!file.good()){
		std::cerr << "File not found: " << filePath << '\n';
		return false;
	}
	return true;
}

/**
 * @file         "FIleRead.cpp"
 * @brief        文件读取时进行跳行操作
 * @param[in]    file
 * @param[in]    linesToSkip
 * @return       bool  true:表示跳过成功  false:跳过失败，一般为到文件尾
 * @note         目前检验比getline();快20%的性能
 * @par History
 *               2024/12/13, Ken Ga  new \n
 */
bool FILEOpt::SkipLines(std::ifstream& file,const int linesToSkip)
{
	 //保存当前位置
	std::streampos currentPos = file.tellg();

	// 跳过linesToSkip行
	int linesSkipped = 0;
	char ch;
	while (linesSkipped < linesToSkip) {
		// 读取下一个字符
		if (!file.get(ch)) {
			return false;  // 如果到达文件末尾，返回false
		}

		// 如果遇到换行符，说明跳过了一行  同时支持WINDOS和UNIX系统
		if (ch == '\n' || ch == '\r') {
			if (ch == '\r' && file.peek() == '\n') {
				file.get();  // 跳过 '\n'
			}
			linesSkipped++;
		}

	}

	// 返回跳过操作是否成功
	return true;
}

std::string FILEOpt::CreatFolder(const std::string& path,const std::string Foldername)
{
	std::filesystem::path posDir = std::filesystem::path(path) / Foldername;
	std::filesystem::create_directories(posDir);
	return posDir.string();
}

void FILEOpt::changeFileExtension(const std::string& filePath, const std::string& newExtension)
{
	// 将文件路径解析为 std::filesystem::path
	std::filesystem::path oldPath(filePath);

	// 检查文件是否存在
	if (!std::filesystem::exists(oldPath)) {
		throw std::runtime_error("File does not exist: " + filePath);
	}

	// 创建新的文件路径
	std::filesystem::path newPath = oldPath;
	newPath.replace_extension(newExtension);

	// 重命名文件
	std::filesystem::rename(oldPath, newPath);

	std::cout << "File renamed from " << oldPath << " to " << newPath << "\n";
}

/**
 * @file         "FIleRead.cpp"
 * @brief        复制文件到指定文件夹
 * @param[in]    fpath
 * @param[in]    oFolder
 * @return       bool
 * @note         特别注意的是，最终的复制的文件名为：源文件夹名@源文件名
 * @par History
 *               2024/12/24, Ken Ga  new \n
 */
bool FILEOpt::copyFileToFolder(const std::string& fpath, const std::string& oFolder)
{
	std::filesystem::path sourceFile(fpath);

	// 检查源文件是否存在
	if (!std::filesystem::exists(sourceFile)) {
		return false; // 源文件不存在，直接返回
	}

	// 提取原文件名和上级文件夹名
	std::string originalFileName = sourceFile.filename().string();
	std::string parentFolderName = sourceFile.parent_path().filename().string();

	// 目标文件夹路径
	std::filesystem::path targetFolder(oFolder);
	if (!std::filesystem::exists(oFolder)) {
		std::filesystem::create_directories(oFolder);
	}

	// 创建新的文件名：<上级文件夹名>@<原文件名>
	std::string newFileName = parentFolderName + "@" + originalFileName;

	// 构建目标文件路径
	std::filesystem::path targetFile = targetFolder / newFileName;

	std::filesystem::copy_file(sourceFile, targetFile, std::filesystem::copy_options::overwrite_existing);
	return true; // 成功复制文件
}

uintmax_t FILEOpt::getFileSize(const std::string& file)
{
	return std::filesystem::file_size(file);
}

std::string FILEOpt::changeStrExtension(const std::string& filePath, const std::string& extension)
{
	size_t dotPosition = filePath.rfind('.');
	if (dotPosition == std::string::npos) {
		// No extension found, append ".pos" to the file path
		return filePath + extension;
	}
	// Replace the extension with ".pos"
	return filePath.substr(0, dotPosition) + extension;
}

void FILEOpt::copyBinFile(const std::string& inputFilePath, const std::string& outputFilePath)
{
	// 打开输入文件，使用二进制模式
	std::ifstream inputFile(inputFilePath, std::ios::binary);
	if (!inputFile.is_open()) {
		throw std::runtime_error("Failed to open input file: " + inputFilePath);
	}

	// 打开输出文件，使用二进制模式
	std::ofstream outputFile(outputFilePath, std::ios::binary);
	if (!outputFile.is_open()) {
		throw std::runtime_error("Failed to open output file: " + outputFilePath);
	}

	// 使用缓冲区分块复制文件
	constexpr size_t bufferSize = 8192; // 8 KB 缓冲区
	char buffer[bufferSize];

	while (inputFile.read(buffer, bufferSize)) {
		// 读取了完整的缓冲区大小
		outputFile.write(buffer, bufferSize);
	}

	// 处理文件的最后一块数据（小于缓冲区大小）
	if (inputFile.gcount() > 0) {
		outputFile.write(buffer, inputFile.gcount());
	}

	// 关闭文件
	inputFile.close();
	outputFile.close();
}




/**
 * @file         "FIleRead.cpp"
 * @brief        递归调用，慎用！！！ 遍历文件夹并对文件夹和文件应用限制条件
 * @param[in]    folderPath
 * @param[in]    fileOperation
 * @return       bool
 * @note         遍历文件夹并对文件夹和文件应用限制条件
 * @par History
 *               2025/01/16, Ken Ga  new \n
 */
bool FILEOpt::traverseALLFolder(const std::string& folderPath,const std::function<void(const std::string&)>& fileOperation) 
{
	std::filesystem::path path(folderPath);
	if (std::filesystem::exists(path) && std::filesystem::is_directory(path)) {
		for (const auto& entry : std::filesystem::directory_iterator(path)) {
			if (std::filesystem::is_directory(entry)) {
				// 文件夹限制：名字包含指定字符串
				if (entry.path().filename().string().find(FolderIdntifier) != std::string::npos) {
					std::cout << "Traversing Directory: " << entry.path() << std::endl;
					Foldercount++;
					traverseALLFolder(entry.path().string(), fileOperation); // 递归调用
				}
			}
			else if (std::filesystem::is_regular_file(entry)) {
				const std::string filename = entry.path().filename().string();
				// 文件限制：扩展名或名字包含指定字符串
				if (filename.find(FileIdentifier) != std::string::npos && entry.path().extension() == Filesuffix) {
					std::cout << "File: " << entry.path() << std::endl;
					Filecount++;
					// 调用传入的函数处理文件
					fileOperation(entry.path().string());  //传入任意void 类型的string参 函数都行
				}
			}
		}
	}
	else {
		std::cerr << "The path is not a valid folder: " << folderPath << std::endl;
		return false;
	}
	return true;
}






//static const char SimbolIdt[] = "$RTK_RAW,";   //$RTK_BASE,  "$RTK_RAW,"
// 将十六进制字符转换为对应的整数值
static unsigned char hexCharToByte(char c) {
	if (c >= '0' && c <= '9') return c - '0';
	if (c >= 'a' && c <= 'f') return 10 + c - 'a';
	if (c >= 'A' && c <= 'F') return 10 + c - 'A';
	throw std::invalid_argument("Invalid hex character");
}

// 将十六进制字符串转换为二进制数据
static std::vector<unsigned char> hexStringToBinary(const std::string& hexStr) {
	std::vector<unsigned char> binaryData;
	for (size_t i = 0; i < hexStr.length(); i += 2) {
		unsigned char highNibble = hexCharToByte(hexStr[i]);
		unsigned char lowNibble = hexCharToByte(hexStr[i + 1]);
		unsigned char byte = (highNibble << 4) | lowNibble;
		binaryData.push_back(byte);
	}
	return binaryData;
}

// 处理 $RTK_RAW 数据并转换为二进制
static std::vector<unsigned char> processRTKRAWData(const std::string& rtkRawData, const char SimbolIdt[]) {
	// 去除前缀 "$RTK_RAW,"
	size_t prefixPos = rtkRawData.find(SimbolIdt);

	int SimbolIdtLen = strlen(SimbolIdt); // 直接使用 strlen;//sizeof(SimbolIdt) / sizeof(char)-1;
	if (prefixPos == std::string::npos) {
		throw std::invalid_argument("Invalid $RTK_RAW data format");
	}
	std::string hexStr = rtkRawData.substr(prefixPos + SimbolIdtLen);

	// 去除校验码（如果有）
	size_t checksumPos = hexStr.find_last_of('*');
	if (checksumPos != std::string::npos) {
		hexStr = hexStr.substr(0, checksumPos);
	}

	// 将十六进制字符串转换为二进制
	return hexStringToBinary(hexStr);
}


void extractRoveFromRAW(const std::string& inputFilePath, const std::string& outputFilePath, const char SimbolIdt[]) {
	std::ifstream inputFile(inputFilePath);
	if (!inputFile.is_open()) {
		throw std::runtime_error("Failed to open input file: " + inputFilePath);
	}

	//std::ofstream outputFile(outputFilePath);
	//if (!outputFile.is_open()) {
	//	throw std::runtime_error("Failed to open output file: " + outputFilePath);
	//}

	std::string line;
	std::vector<unsigned char> binaryData;
	std::ofstream outputFile(outputFilePath, std::ios::binary);
	if (!outputFile.is_open())printf("%s", "输出文件无法打开");
	while (std::getline(inputFile, line)) {
		// Check if the line starts with "$GNGGA"
		//if (line.find(SimbolIdt) == 0) {
		if (line.find(SimbolIdt) != std::string::npos) {
			//size_t checksumPos = line.find_last_of('*');
			//std::cout << line.substr(9, checksumPos-3) << std::endl;
			//outputFile << line.substr(9, checksumPos);// << '\n';
			binaryData=processRTKRAWData(line,SimbolIdt);
			outputFile.write(reinterpret_cast<const char*>(binaryData.data()), binaryData.size());
		}
	}
	outputFile.close();

	//constexpr size_t bufferSize = 4096; // 4 KB 缓冲区
	//char buffer[bufferSize];

	//while (inputFile.read(buffer, bufferSize)) {
	//	outputFile.write(buffer, inputFile.gcount());
	//}

	//// 处理最后一块（如果有剩余）
	//outputFile.write(buffer, inputFile.gcount());


	// Close files
	inputFile.close();
	//outputFile.close();
}

/**
 * @file         "FIleRead.cpp"
 * @brief        将文件改名，并移动到上一级目录
 * @param[in]    oldFilePath
 * @return       bool
 * @note         
 * @par History
 *               2024/12/27, Ken Ga  new \n
 */
bool renameTHObsFile(const std::string& oldFilePath) {
	// 提取文件名和文件后缀
	std::string fileName = std::filesystem::path(oldFilePath).filename().string();
	std::string fileExtension = std::filesystem::path(oldFilePath).extension().string();

	// 确保文件扩展名是 .24o
	if (!(fileExtension == ".25o" || fileExtension==".25p")) {
		return false;
	}

	// 文件名格式：r20241206001.24o，提取日期和时段
	//if (fileName.size() < 15 || fileName[0] != 'r') {
	//	std::cerr << "Invalid file name format." << std::endl;
	//	return false;
	//}

	std::string date = fileName.substr(0, 8);  // 20241206
	std::string timePeriod = fileName.substr(10, 2);  // 01 或 02

	// 根据日期转换为目标格式（假设这里是手动调整）
	// 比如 20241206 -> 20240913
	std::string newDate = date.substr(0, 4) + date.substr(4, 2) + date.substr(6, 2);  // 假设没有更复杂的日期转换需求

	// 创建新文件名
	std::string newFileName = "ROVE_" + newDate + "_"+timePeriod + "_TH1100_" + "01" + fileExtension;

	//// 获取新的文件路径
	//std::filesystem::path newFilePath = std::filesystem::path(oldFilePath).parent_path() / newFileName;
	// 获取文件的当前路径
	std::filesystem::path currentDir = std::filesystem::path(oldFilePath).parent_path();
	// 获取上一级目录路径
	std::filesystem::path parentDir = currentDir.parent_path();
	// 获取新的文件路径（在上一级目录中）
	std::filesystem::path newFilePath = parentDir / newFileName;

	try {
		// 重命名文件
		std::filesystem::rename(oldFilePath, newFilePath);
		return true;
	}
	catch (const std::exception& e) {
		return false;
	}
}



// 结构体用来存储文件的日期、时间、以及原始路径
struct FileInfo {
	std::string date;  // 格式: YYYYMMDD
	std::string time;  // 格式: HHMMSS
	std::string filePath;
	std::string extision;

	// 用于排序：先按日期，再按时间排序
	bool operator<(const FileInfo& other) const {
		if (extision !=other.extision) {
			if (extision == ".24o")return true;
			else return false;
		}
		if (date != other.date)
			return date < other.date;
		return time < other.time;
	}
};

// 将时间从UTC+8转换为UTC
static std::string convertToUTC(const std::string& date, const std::string& time) {
	// 解析日期和时间
	int year = std::stoi(date.substr(0, 4));
	int month = std::stoi(date.substr(4, 2));
	int day = std::stoi(date.substr(6, 2));
	int hour = std::stoi(time.substr(0, 2));
	int minute = std::stoi(time.substr(2, 2));
	int second = std::stoi(time.substr(4, 2));

	// 将日期和时间存入 std::tm
	std::tm tm = {};
	tm.tm_year = year - 1900; // tm.tm_year 是从 1900 年开始的
	tm.tm_mon = month - 1;    // tm.tm_mon 是从 0 开始的
	tm.tm_mday = day;
	tm.tm_hour = hour;
	tm.tm_min = minute;
	tm.tm_sec = second;

	// 将时间从 UTC+8 转换为 UTC（减去 8 小时）
	tm.tm_hour -= 0;  //8

	// 如果小时数变为负数，则需要调整日期
	if (tm.tm_hour < 0) {
		tm.tm_hour += 24;  // 恢复到正值
		tm.tm_mday -= 1;   // 减去一天
		if (tm.tm_mday < 1) {
			// 处理月和年的变化
			tm.tm_mon -= 1;
			if (tm.tm_mon < 0) {
				tm.tm_mon = 11;  // 12月
				tm.tm_year -= 1; // 年份减1
			}

			// 根据月份和年份调整日期
			// 获取该月的天数
			const int daysInMonth[] = { 31, (tm.tm_year % 4 == 0 && (tm.tm_year % 100 != 0 || tm.tm_year % 400 == 0)) ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
			tm.tm_mday = daysInMonth[tm.tm_mon]; // 获取该月的最后一天
		}
	}

	// 格式化为 "YYYYMMDD" 和 "HHMMSS"
	std::ostringstream newDateStream;
	newDateStream << std::setw(4) << std::setfill('0') << tm.tm_year + 1900
		<< std::setw(2) << std::setfill('0') << tm.tm_mon + 1
		<< std::setw(2) << std::setfill('0') << tm.tm_mday;

	std::ostringstream newTimeStream;
	newTimeStream << std::setw(2) << std::setfill('0') << tm.tm_hour
		<< std::setw(2) << std::setfill('0') << tm.tm_min
		<< std::setw(2) << std::setfill('0') << tm.tm_sec;

	return newDateStream.str() + "_" + newTimeStream.str();
}


// 提取日期和时间并构造文件信息
static FileInfo extractFileInfo(const std::string& filePath) {
	FileInfo info;
	std::filesystem::path p(filePath);
	std::string fileName = p.filename().string();

	int date = std::stoi(fileName);
	// 假设文件名格式是：YYYYMMDD_HHMMSS.24o
	info.date = fileName.substr(0, 8);  // 提取日期部分：20241206
	info.time = fileName.substr(9, 6);  // 提取时间部分：111518
	info.filePath = filePath;

	std::string newDateAndTime = convertToUTC(info.date, info.time);
	info.date = newDateAndTime.substr(0, 8);  // 更新为 UTC 日期
	info.time = newDateAndTime.substr(9, 6);  // 更新为 UTC 时间

	return info;
}

// 获取文件夹下所有符合条件的文件路径
std::vector<FileInfo> traverseTHPdir(const std::string& folderPath) {
	std::vector<FileInfo> allFiles;
	for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			const std::filesystem::path& filePath = entry.path();
			//if (filePath.filename().string().size() != 19)continue;  //这里的文件大小需要加上后缀名
			if (filePath.extension() == ".25o"||filePath.extension()==".25p") {  // 只处理 .24o 文件
				if (filePath.filename().string().find("BASE") != std::string::npos)continue;  //基站文件跳过
				FileInfo fileInfo = extractFileInfo(filePath.string());
				fileInfo.extision = filePath.extension().string();
				allFiles.push_back(fileInfo);
			}
		}
	}
	return allFiles;
}

// 用于为同一天的文件分配时间段号
namespace fs = std::filesystem;
static void renameTHObsFiles2(const std::vector<FileInfo>& sortedFiles, const std::string& folderPath,const std::string & fuix) {
	std::string currentDate;
	int timeSlotCounter = 1;  // 用于分配时间段号 01, 02, ...

	for (const auto& fileInfo : sortedFiles) {

		if (fileInfo.extision != fuix)continue;
		// 判断日期是否变化，若变化，重置时间段计数器
		if (fileInfo.date != currentDate) {
			currentDate = fileInfo.date;
			timeSlotCounter = 1;  // 重置时间段计数器
		}
		else {
			++timeSlotCounter;  // 增加时间段号
		}
		
		std::string newFileName;
		std::string newDate = fileInfo.date;
		std::string timeSlot = (timeSlotCounter < 10) ? "0" + std::to_string(timeSlotCounter) : std::to_string(timeSlotCounter);

		// 生成新文件名，格式：ROVE_YYYYMMDD_Tx_TH1200_yx_24o
		newFileName = "ROVE_" + newDate + "_" + timeSlot + "_TH1100_" + "01" + fuix;

		//// 获取原始文件的路径，并生成新文件路径
		//std::filesystem::path newFilePath = fs::path(folderPath) / newFileName;

		// 获取文件的当前路径
		std::filesystem::path currentDir = std::filesystem::path(folderPath);
		// 获取上一级目录路径
		std::filesystem::path parentDir = currentDir.parent_path();
		// 获取新的文件路径（在上一级目录中）
		std::filesystem::path newFilePath = parentDir / newFileName;

		//// 获取文件的当前路径
		//std::filesystem::path currentDir = std::filesystem::path(folderPath);
		//// 获取新的文件路径（在上一级目录中）
		//std::filesystem::path newFilePath = currentDir / newFileName;

		try {
			// 重命名文件
			fs::rename(fileInfo.filePath, newFilePath);
			std::cout << "Renamed: " << fileInfo.filePath << " -> " << newFilePath << std::endl;
		}
		catch (const std::exception& e) {
			std::cerr << "Error renaming file: " << e.what() << std::endl;
		}

		
	}
}

bool PDirTHobs(std::string folderPath)
{

	// 获取文件夹下所有符合条件的文件
	std::vector<FileInfo> files = traverseTHPdir(folderPath);

	// 按日期和时间排序
	std::sort(files.begin(), files.end());

	// 重命名并移动文件
	renameTHObsFiles2(files, folderPath, ".25o");
	renameTHObsFiles2(files, folderPath, ".25p");


	return true;
}



void renameTH_OneWallsDATA(std::string &filepath)
{
	// 解析文件名

	std::filesystem::path path = filepath;
	std::string oldFileName = path.filename().string();  // 获取文件名

	// 查找并提取文件名中的日期部分和编号部分
	size_t startPos = oldFileName.find('r');  // 假设文件名以'r'开始

	if (startPos != std::string::npos) {
		std::string datePart = oldFileName.substr(startPos + 1, 8);  // 日期部分（20241214）
		std::string numPart = oldFileName.substr(startPos + 9, 3);   // 编号部分（007）

		// 构建新的文件名
		std::string newFileName = datePart + "_" + numPart + ".rtcm3";

		// 构造新路径
		fs::path newFilePath = path.parent_path() / newFileName;

		// 输出修改前后的路径
		std::cout << "Renaming: " << filepath << " -> " << newFilePath.string() << std::endl;

		// 使用 std::filesystem 重命名文件
		std::filesystem::rename(filepath, newFilePath.string());
			
	}
	else {
		std::cerr << "Invalid filename format!" << std::endl;
	}
}


namespace fs = std::filesystem;

// 提取文件日期，例如 "ROVE_20241203_01_TH1200_01.24o" -> "20241203"
std::string extractDate(const std::string& fileName) {
	// 匹配文件名中的日期部分，例如 "BASE_202412050516_TH_01.24o" -> "20241205"
	std::regex dateRegex(R"(_(\d{8}))");
	std::smatch match;
	if (std::regex_search(fileName, match, dateRegex)) {
		return match[1]; // 返回匹配的日期部分
	}
	return ""; // 如果没有匹配到，则返回空字符串
}


// 提取基准站文件名，与流动站日期匹配
static std::string getBaseFile(const std::string& roveDate, const std::map<std::string, std::string>& baseFiles) {
	auto it = baseFiles.find(roveDate);
	if (it != baseFiles.end()) {
		return it->second; // 返回匹配的 BASE 文件名
	}
	return "";
}

// 生成 .opt 文件
static void generateOptFile(const std::string& roveFile, const std::string& baseFile, const std::string& outputPath) {
	std::string roveFileName = fs::path(roveFile).filename().string();
	std::string roveFilePath = fs::path(roveFile).parent_path().string();

	// 生成星历文件名，将 .24o 后缀改为 .24p
	std::string naviFileName = roveFileName;
	naviFileName.replace(roveFileName.find(".25o"), 4, ".25p");

	// 生成 opt 文件路径
	std::string optFileName = roveFileName;
	optFileName.replace(optFileName.find(".25o"), 4, ".opt");
	std::string optFilePath = outputPath + "/" + optFileName;

	// 定义要替换的三行内容
	std::string baseLine = "GNSS_BaseFile =                  " + baseFile + "                  ; <base station fileName>";
	std::string roveLine = "GNSS_RoveFile =                  " + roveFileName + "                 ; <rove station fileName>";
	std::string naviLine = "GNSS_NaviFile =                  " + naviFileName + "                 ; <Eph filename>";

	// 打开模板文件 (假设模板文件为 "template.opt"，在同一目录下)
	std::ifstream templateFile("template.opt");
	if (!templateFile) {
		std::cerr << "无法打开模板文件 template.opt，请确保模板文件存在！" << std::endl;
		return;
	}

	// 读取模板文件内容
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(templateFile, line)) {
		lines.push_back(line);
	}
	templateFile.close();

	// 修改模板文件中的对应行
	for (auto& l : lines) {
		if (l.find("GNSS_BaseFile = ") != std::string::npos) {
			l = baseLine;
		}
		else if (l.find("GNSS_RoveFile = ") != std::string::npos) {
			l = roveLine;
		}
		else if (l.find("GNSS_NaviFile = ") != std::string::npos) {
			l = naviLine;
		}
	}

	// 将修改后的内容写入新的 opt 文件
	std::ofstream optFile(optFilePath);
	if (!optFile) {
		std::cerr << "无法创建文件: " << optFilePath << std::endl;
		return;
	}

	for (const auto& l : lines) {
		optFile << l << "\n";
	}
	optFile.close();

	std::cout << "生成文件: " << optFilePath << std::endl;
}


void mergeBASEObsFiles(const std::filesystem::path& dir) {
	// 匹配模式：BASE_日期_时段_设备型号.250
	const std::regex reg(R"(BASE_(\d{8})_(\d{2})_([^.]+)\.25o)");
	std::map<std::string, std::vector<fs::path>> groups; // 日期+设备作为key

	// 单次遍历完成文件收集和分组
	for (const auto& entry : fs::directory_iterator(dir)) {
		std::cmatch match;
		const std::string name = entry.path().filename().string();

		if (std::regex_match(name.c_str(), match, reg)) {
			groups[match[1].str() + "_" + match[3].str()].push_back(entry.path());
		}
	}

	// 处理每个分组
	for (const auto& [key, files] : groups) {
		// 按时段排序
		std::vector<fs::path> sorted(files);
		std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) {
			return a.filename().string() < b.filename().string();
		});

		// 创建输出文件
		const fs::path out = dir / ("BASE_" + key + ".25o");
		std::ofstream fout(out);
		if (!fout) {
			std::cerr << "创建文件失败: " << out << std::endl;
			continue;
		}

		// 合并内容
		bool isFirst = true;
		for (const auto& file : sorted) {
			std::ifstream fin(file);
			if (!fin) {
				std::cerr << "打开失败: " << file << std::endl;
				continue;
			}

			if (isFirst) {
				fout << fin.rdbuf();
				isFirst = false;
			}
			else {
				//skipRinexHead(fin); // 使用现有的跳过头部函数
				fout << fin.rdbuf();
			}
			std::cout << "已合并: " << file.filename() << std::endl;
		}

		std::cout << "生成合并文件: " << out.filename() << "\n\n";
	}
}


// 主逻辑函数
void setTHOptfile(const std::string& folderPath) {
	// 保存所有 BASE 和 ROVE 文件
	std::map<std::string, std::string> baseFiles; // key 为日期 YYYYMMDD
	std::vector<std::string> roveFiles;

	// 遍历文件夹
	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			std::string fileName = entry.path().filename().string();
			if (fileName.find("BASE") != std::string::npos && fileName.find(".25o") != std::string::npos) {
				// 提取日期并存储 BASE 文件
				std::string date = extractDate(fileName);
				if (!date.empty()) {
					baseFiles[date] = fileName;
				}
			}
			else if (fileName.find("ROVE") != std::string::npos && fileName.find(".25o") != std::string::npos) {
				// 存储 ROVE 文件
				roveFiles.push_back(entry.path().string());
			}
		}
	}

	// 输出路径（重新设置opt文件夹）
	FILEOpt::CreatFolder(folderPath, "opt");
	std::string outputPath = folderPath+"\\"+"opt";

	// 为每个 ROVE 文件生成 .opt 文件
	for (const auto& roveFile : roveFiles) {
		std::string roveDate = extractDate(roveFile);
		std::string baseFile = getBaseFile(roveDate, baseFiles);
		if (!baseFile.empty()) {
			generateOptFile(roveFile, baseFile, outputPath);
		}
		else {
			std::cerr << "未找到匹配的 BASE 文件: " << roveFile << std::endl;
		}
	}
}




// 提取 .24o 文件中的时间信息并生成新的文件名
static std::string extractTimeFromRinexObs(const fs::path& filePath) {
	std::ifstream file(filePath);
	if (!file.is_open()) {
		std::cerr << "无法打开文件: " << filePath << std::endl;
		return "";
	}

	std::string line;
	bool headerEnded = false;
	std::regex timePattern(R"(^>\s*(\d{4})\s+(\d{2})\s+(\d{2})\s+(\d{2})\s+(\d{2})\s+([\d.]+))");

	while (std::getline(file, line)) {
		// 检测 END OF HEADER
		if (!headerEnded && line.find("END OF HEADER") != std::string::npos) {
			headerEnded = true;
			continue;
		}

		// 如果已到 HEADER 后，匹配时间行
		if (headerEnded && line[0] == '>') {
			std::smatch match;
			if (std::regex_search(line, match, timePattern)) {
				std::ostringstream newName;
				newName << match[1] << match[2] << match[3] << "_"
					<< match[4] << match[5] << match[6].str().substr(0, 2); // 精确到秒
				return newName.str();
			}
		}
	}

	std::cerr << "未能在文件中找到有效的时间信息: " << filePath << std::endl;
	return "";
}

// 处理目标文件夹，重命名ROVE .24o 和 .24p 文件
void RenameObsEphFILE_ROVE(const std::string& folderPath) {
	fs::path targetDir(folderPath);

	if (!fs::exists(targetDir) || !fs::is_directory(targetDir)) {
		std::cerr << "无效的目录路径: " << folderPath << std::endl;
		return;
	}

	for (const auto& entry : fs::directory_iterator(targetDir)) {
		if (entry.is_regular_file() && entry.path().extension() == ".25o" ) {
			if (entry.path().filename().string().find("BASE") != std::string::npos)
				continue; // 检测到 BASE 文件夹跳过
			fs::path oFile = entry.path();
			std::string baseName = extractTimeFromRinexObs(oFile);

			if (!baseName.empty()) {
				// 新的 .24o 文件名
				fs::path newOFile = oFile.parent_path() / (baseName + ".25o");

				try {
					// 重命名 .24o 文件
					fs::rename(oFile, newOFile);
					std::cout << "已重命名: " << oFile.filename() << " -> " << newOFile.filename() << std::endl;

					// 找到对应的 .24p 文件并重命名
					fs::path pFile = oFile;
					pFile.replace_extension(".25p");
					if (fs::exists(pFile)) {
						fs::path newPFile = pFile.parent_path() / (baseName + ".25p");
						fs::rename(pFile, newPFile);
						std::cout << "已重命名: " << pFile.filename() << " -> " << newPFile.filename() << std::endl;
					}
				}
				catch (const fs::filesystem_error& e) {
					std::cerr << "重命名文件时发生错误: " << e.what() << std::endl;
				}
			}
		}
	}
}

//static const double gs_WGS84_a = 6378137.0;                         ///< earth semimajor axis (WGS84) (m)
static const double gs_WGS84_e2 = 0.00669437999014;



void LLH2XYZ(const double LLH[3], double XYZ[3])
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;


	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	double v = a / sqrt(1.0 - e2 * sinp * sinp);

	XYZ[0] = (v + LLH[2]) * cosp * cosl;
	XYZ[1] = (v + LLH[2]) * cosp * sinl;
	XYZ[2] = (v * (1.0 - e2) + LLH[2]) * sinp;
}


double DegMinSec2D(const double& X){
	int Deg = int(X);
	double min = ((X - Deg)*100);
	return Deg + min / 60;
}


void XYZ2LLH(const double XYZ[3], double LLH[3])
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	double X = XYZ[0];
	double Y = XYZ[1];
	double Z = XYZ[2];

	double r2 = X * X + Y * Y;
	double z = 0.0;
	double zk = 0.0;
	double v = a;
	double sinp = 0.0;

	for (z = Z, zk = 0.0; fabs(z - zk) >= 1E-4;)
	{
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = a / sqrt(1.0 - e2 * sinp * sinp);
		z = Z + v * e2 * sinp;
	}
	LLH[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (Z > 0.0 ? PI / 2.0 : -PI / 2.0);
	LLH[1] = r2 > 1E-12 ? atan2(Y, X) : 0.0;
	LLH[2] = sqrt(r2 + z * z) - v;
}




#include <numeric>  // for accumulate
#include <algorithm>
using namespace std;
// 计算均值
static double mean(const vector<double>& data) {
	if (data.empty()) return 0.0;
	return accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// 计算标准差
static double standardDeviation(const vector<double>& data, double meanValue) {
	if (data.size() < 2) return 0.0;

	double sumSq = 0.0;
	for (double val : data) {
		sumSq += (val - meanValue) * (val - meanValue);
	}
	return sqrt(sumSq / (data.size() - 1)); // 样本标准差 (n-1)
}

// 计算去除3sigma后的均值
double meanWithoutOutliers(const vector<double>& data) {
	if (data.size() < 3) return mean(data); // 数据太少，直接返回原均值

	double mu = mean(data);
	double sigma = standardDeviation(data, mu);

	// 过滤掉超出 mean ± 3*sigma 的点
	vector<double> filteredData;
	for (double val : data) {
		if (fabs(val - mu) <= 3 * sigma) {
			filteredData.push_back(val);
		}
	}

	return mean(filteredData); // 计算新均值
}



void getRefPos(int PosType){
	std::string folderPath;
	std::cout << "请输入文件路径:定位的pos文件夹路径)" << std::endl << std::endl;
	std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径

	std::ifstream inputFile(folderPath);
	if (!inputFile.is_open()) {
		throw std::runtime_error("Failed to open input file: " + folderPath);
	}

	std::string line;
	std::string token;
	double XYZ[3]{}, BLH[3]{};
	int solStatus = 0;
	std::vector<double>Pos[3];
	double Rad = PI / 180;
	while (std::getline(inputFile, line)) {
		// Check if the line starts with "$GNGGA"
		if (line.find("$GNGGA") == 0 && PosType==0) {

			std::stringstream ss(line);
			std::vector<std::string> tokens;
			while (std::getline(ss, token, ',')) {
				tokens.push_back(token);
			}
			if (tokens[6] != "4")continue;

			BLH[0] = std::stod(tokens[2]) / 100.0;
			BLH[1] = std::stod(tokens[4]) / 100.0;
			BLH[2] = std::stod(tokens[9]);
			for (size_t i = 0; i < 2; i++)BLH[i] = DegMinSec2D(BLH[i]) * Rad;  //度分转为度

			if (tokens[3] == "S")BLH[0] = -BLH[0];
			if (tokens[5] == "W")BLH[1] = -BLH[1];
			LLH2XYZ(BLH, XYZ);
			for (size_t i = 0; i < 3; i++)Pos[i].push_back(XYZ[i]);
		}
		// IE 的HW格式
		if (PosType == 1 && line.find("Fixed") != std::string::npos) {
			std::stringstream ss(line);
			std::vector<std::string> tokens;
			tokens = FILEOpt::splitString(line, ' ');
			if (tokens[5] != "1")continue;
			for (size_t i = 0; i < 3; i++)XYZ[i] = stod(tokens[16 + i]);
			for (size_t i = 0; i < 3; i++)Pos[i].push_back(XYZ[i]);
		}
	}
	if (Pos[0].size() == 0)return ;

	double RefPos[3]{};
	for (int j = 0; j < 3; j++)RefPos[j] = meanWithoutOutliers(Pos[j]);
	//for (size_t i = 0; i < Pos[0].size(); i++)
	//{
	//    for (int j = 0; j < 3; j++)RefPos[j] += Pos[j][i];
	//}
	//for (int j = 0; j < 3; j++)RefPos[j] = RefPos[j]/ Pos[0].size();

	std::cout << "该文件固定解的平均值为 XYZ： " << std::endl;
	printf("%12.4f,%12.4f,%12.4f\n", RefPos[0], RefPos[1], RefPos[2]);
	XYZ2LLH(RefPos, BLH);

	std::cout << "该文件固定解的平均值为 BLH： " << std::endl;
	printf("%12.9f,%13.9f,%8.4f\n", BLH[0] * 180 / PI, BLH[1] * 180 / PI, BLH[2]);

	printf("\n\n");
}



const double Rad = PI / 180;
static bool getIPSpos(const std::string& line, double IPSdata[],UTCTIME &utc,GPSTIME &gpsIPS) {
	double BLH[3], XYZ[3];
	std::stringstream ss(line);
	std::string token;
	std::vector<std::string> tokens;
	bool status;
	while (std::getline(ss, token, ',')) {
		tokens.push_back(token);
	}

	static bool isCorssDay = false;
	utc.Hour = std::stoi(tokens[1].substr(0, 2));
	utc.Minute = std::stoi(tokens[1].substr(2, 2));
	utc.Second = std::stod(tokens[1].substr(4));
	if (utc.Hour == 0 && utc.Minute == 0 && utc.Second == 0.0&& !isCorssDay)
	{
		utc.Day += 1;
		isCorssDay = true;
	}
	BLH[0] = std::stod(tokens[2]) / 100.0;
	BLH[1] = std::stod(tokens[4]) / 100.0;
	BLH[2] = std::stod(tokens[9]);
	for (size_t i = 0; i < 2; i++)BLH[i] = DegMinSec2D(BLH[i]) * Rad;  //度分转为度
	gpsIPS = UTCtoGPST(utc);
	if (tokens[3] == "S")BLH[0] = -BLH[0];
	if (tokens[5] == "W")BLH[1] = -BLH[1];
	LLH2XYZ(BLH, XYZ);
	IPSdata[0] = gpsIPS.getSow();
	if (gpsIPS.secWeek >= 259218) {

		int a = 0;
	}
	for (size_t i = 0; i < 3; i++)
	{
		//IPSdata[i + 1] = XYZ[i];
		IPSdata[i + 1] = BLH[i];
	}
	std::stoi(tokens[6]) == 4 ? status = true : status = false;
	return status;
}


static bool getRefStationChangeTime() {
	std::string folderPath;
	std::cout << "请输入文件路径:    (原始观测数据的文件夹路径，如POSRTK.log)" << std::endl;
	std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径
	
}


// 输出ecefXYZ 和姿态(deg)
static bool extraDeSayRealPos(const std::string& line, double Data[]) {
	double BLH[3]{};
	vector<std::string>raw = FILEOpt::splitString(line, ',');
	Data[0] = std::stod(raw[6]);  // 时间(周内秒)

	BLH[0] = std::stod(raw[11]);	// 位置
	BLH[1] = std::stod(raw[12]);
	BLH[2] = std::stod(raw[13]);

	//CoordinateTransform::BLH2XYZ(BLH, Data + 1);
	for (size_t i = 0; i < 3; i++)Data[i + 1] = BLH[i];

	Data[6] = std::stod(raw[18]);  //姿态  航向、俯仰、横滚
	Data[5] = std::stod(raw[19]);
	Data[4] = std::stod(raw[20]);

	return 1;
}


static bool extraDaoYuanPos(const std::string& line, double Data[]) {
	std::stringstream ss(line);
	std::string token;
	std::vector<std::string> tokens;
	double XYZ[3];
	while (std::getline(ss, token, ',')) {
		tokens.push_back(token);
	}
	Data[0] = std::stod(tokens[0]);
	Data[1] = std::stod(tokens[2]);//* Rad;  //纬度
	Data[2] = std::stod(tokens[1]);//*Rad;  //经度
	Data[3] = std::stod(tokens[3]);

	Data[4] = std::stod(tokens[7]);  //航向，俯仰，横滚   输入：横滚，俯仰，航向
	Data[5] = std::stod(tokens[8]);  // 俯仰
	Data[6] = std::stod(tokens[9]);  // 横滚
	//LLH2XYZ(Data + 1, Data + 1);
	return true;
}


static bool extraIEPos(const std::string& line, double Data[],int skipline=24) {
	std::stringstream ss(line);
	std::string token;
	std::vector<std::string> tokens;
	double XYZ[3];
	tokens = FILEOpt::splitString(line, ' ');
	Data[0] = std::stod(tokens[1]);

	//Data[1] = std::stod(tokens[2])* Rad;  //纬度
	//Data[2] = std::stod(tokens[3])*Rad;  //经度
	//Data[3] = std::stod(tokens[4]);     // 高

	Data[1] = std::stod(tokens[16]);  //X
	Data[2] = std::stod(tokens[17]);  //Y
	Data[3] = std::stod(tokens[18]);


	//Data[4] = std::stod(tokens[7]);  //航向，俯仰，横滚   输入：横滚，俯仰，航向
	//Data[5] = std::stod(tokens[8]);  // 俯仰
	//Data[6] = std::stod(tokens[9]);  // 横滚

	if (tokens[6] == "Fixed")return true;
	else return false;
	//LLH2XYZ(Data + 1, Data + 1);
}

double Distance2D(const double &a, const double &b) {
	return sqrt(pow(a, 2) + pow(b, 2));
}


/**
 * @file         "FIleRead.cpp"
 * @brief        
 * @return       void
 * @note         注意需要修改对应的时间！！！
 * @par History
 *               2025/03/7, Ken Ga  new \n
 */
void DifPos(/*std::string &InputIPS,std::string &InputCouple*/) {
	const double XYZref[] = { -2162417.3932,4394659.4002,4072012.9125 };

	double enuIPS[3]{}, enuCoup[3]{};
	//std::string InputIPS = "E:\\CDTH\\2025@Couple\\20250516\\20250516_100C_POS\\516_605_100C_TX.asc";  "E:\CDTH\OtherVision\POSRTK_Car\PredictVelPOS.pos"
	std::string InputIPS = "F:\\CDTH\\202406@CDTH@CropperGNSS\\KG\\20241203@TH1100@BSDVRS@UGV@Trees\\res\\POSRTKDLLV34Tmp_ROVE_20241202_19_TH1200_01\\SNR@EL@Var@COM@POSRTK.pos";
	//std::string InputIPS = "E:\\CDTH\\2025@Couple\\20250516\\res\\POSRTKDLLV337_Rove\\POSRTK_L.pos";
	std::string InputCouple = "F:\\CDTH\\202406@CDTH@CropperGNSS\\KG\\20241203@TH1100@BSDVRS@UGV@Trees\\res\\POSRTKDLLV34Tmp_ROVE_20241202_19_TH1200_01\\FLO_POSRTK.pos";  //IE或者DaoYuan组合导航
	//std::string InputCouple = "E:\\CDTH\\2025@Couple\\20250516\\res\\POSRTKDLLV336QF_Rove\\IE.pos";
	std::string OutputfilePath = FILEOpt::getParentDir(InputIPS) + "\\FLO" + ".sta";
	std::string OutputfileNorm = FILEOpt::getParentDir(InputIPS) +"\\FLO" + ".txt";

	std::ifstream iFileIPS(InputIPS);
	std::ifstream iFileCoup(InputCouple);
	if (!iFileIPS.is_open() || !iFileCoup.is_open()) {
		throw std::runtime_error("Failed to open input file"); 
	}
	UTCTIME utc; GPSTIME gpsIPS, gpsCoup;
	std::string line;
	std::string token;
	double XYZ[3]{}, BLH[3]{}, newBLH[3]{};
	int ahead = 0;
	
	// 后面三个为杆臂
	double IPSdata[10]{ 0.0,0.0,0.0,0.0,0.0,0.0,0.0 ,349e-3,-1251e-3,-1478e-3 };
	double  Coupdata[10] = { 0.0,-2162417.3932,4394659.4002,4072012.9125,0.0,0.0,0.0 ,-39e-3,-178e-3, -94e-3 };   //  1488e-3,403e-3, -1390e-3
	utc.Year = 2024; utc.Month = 12; utc.Day = 2;
	
	char status = 0x00;
	//354083
	while (1) {
		if (ahead == 0 || ahead == 1) {		    //	0代表 时间同步，1代表IPS结果落后
			if (!std::getline(iFileIPS, line))break;
			//Check if the line starts with "$GNGGA"
			if (line.find("$GNGGA") != 0)continue;
			getIPSpos(line, IPSdata, utc, gpsIPS);

			//std::vector<string>data = FILEOpt::splitString(line, ',');
			//for (size_t i = 0; i < 4; i++)IPSdata[i] = std::stod(data[i]);

			//static int skiplineIPS = 14;  // DaoYuan 26
			//if (skiplineIPS > 0) {
			//	std::getline(iFileIPS, line);
			//	skiplineIPS--;
			//	continue;
			//}
			//std::vector<string>data = FILEOpt::splitString(line, ' ');
			//for (size_t i = 0; i < 4; i++)IPSdata[i] = std::stod(data[i+1]);


			////status=extraDeSayRealPos(line, IPSdata);
			//// 杆臂及坐标转化
			//if (0) {
			//	ins_to_gnss(IPSdata + 1, IPSdata + 4, IPSdata + 7, newBLH);
			//	CoordinateTransform::BLH2XYZ(newBLH, IPSdata + 1, false);
			//}
		}


		if (ahead == 0 || ahead == -1) {		//	0代表 时间同步，-1代表组合导航结果落后

			//if (!std::getline(iFileCoup, line))break;
			//static int skiplineCouple = 24;  // DaoYuan 26

			//if (skiplineCouple > 0) {
			//	skiplineCouple--;
			//	ahead = -1;
			//	continue;
			//}
			//if (line == "")break;

			////extraIEPos(line, Coupdata);
			//extraDaoYuanPos(line, Coupdata);
			////extraDeSayRealPos(line, Coupdata);
			//Coupdata[0] = (int)(Coupdata[0] * 1e3)* 1.0*1e-3;

			//// 杆臂及坐标转化
			//if (1) {
			//	ins_to_gnss(Coupdata + 1, Coupdata + 4, Coupdata + 7, newBLH);
			//	CoordinateTransform::BLH2XYZ(newBLH, Coupdata + 1,false);
			//}

			if (!std::getline(iFileCoup, line))break;
			if (line.find("$GNGGA") != 0)continue;
			status = getIPSpos(line, Coupdata, utc, gpsCoup);

			//std::vector<string>data = FILEOpt::splitString(line, ',');
			//for (size_t i = 0; i < 4; i++)Coupdata[i] = std::stod(data[i]);

			

		}

		//if (Coupdata[0] - IPSdata[0] > 0.002) {  //认为不同步	
		//	ahead = 1;
		//	continue;
		//}

		//else if (IPSdata[0] - Coupdata[0] > 0.002) {  //表示组合导航结果落后
		//	ahead = -1;
		//	continue;
		//}
		//else ahead = 0;
		if (fabs(IPSdata[0] - 123633.0) < 0.001) {
			printf("hit\n");
		}
		//if (status == 0)continue;
		double bIPS[3], bCoup[3];
		CoordinateTransform::Vxyz2enu(XYZref,IPSdata+1, enuIPS,false);
		//TransformNtoB(enuIPS, Coupdata + 4, bIPS);
		CoordinateTransform::Vxyz2enu(XYZref, Coupdata+1, enuCoup,false);
		//TransformNtoB(enuCoup, Coupdata + 4, bCoup);
		static FILE* fpSta = NULL;
		static FILE* fpnorm = NULL;
		if (fpSta == NULL) {
			fopen_s(&fpSta, OutputfilePath.c_str(), "w");
			fopen_s(&fpnorm, OutputfileNorm.c_str(), "w");
		}
		if (fpSta &&ahead == 0) {   // 时间同步了
			//fprintf(fpSta, "%9.2f,%13.10f,%13.10f,%6.3f,%8.5f,%8.5f,%8.5f\n", IPSdata[0], IPSdata[1], IPSdata[2], IPSdata[3],Coupdata[4],Coupdata[5],Coupdata[6]);
			//fprintf(fpSta, "%9.2f,%13.10f,%13.10f,%6.3f\n", Coupdata[0], Coupdata[1], Coupdata[2], Coupdata[3]);
			//fprintf(fpSta, "%8.1f,%5.3f,%5.3f,%5.3f,%d\n", IPSdata[0], (IPSdata[1] - Coupdata[1])* lat_CONSBLH2Meter, (IPSdata[2] - Coupdata[2])* lon_CONSBLH2Meter, IPSdata[3] - Coupdata[3], status);
			fprintf(fpSta, "%8.1f,%5.3f,%5.3f,%5.3f,%d\n", IPSdata[0], enuIPS[0] - enuCoup[0], enuIPS[1] - enuCoup[1], enuIPS[2] - enuCoup[2], status);
			//fprintf(fpSta, "%8.1f,%5.3f,%5.3f,%5.3f,%d,%6.3f\n", IPSdata[0], IPSdata[1] - Coupdata[1], IPSdata[2] - Coupdata[2], IPSdata[3] - Coupdata[3], status, Distance3D(IPSdata + 1, Coupdata + 1));
			printf("%8.1f,%6.3f\n", IPSdata[0], Distance3D(IPSdata + 1, Coupdata + 1));
			//fprintf(fpSta, "%8.1f,%5.3f,%5.3f,%5.3f,%d\n", IPSdata[0], bIPS[0] - bCoup[0], bIPS[1] - bCoup[1], bIPS[2] - bCoup[2], status);
			//printf("%8.1f,%6.3f\n", IPSdata[0], Distance3D(bIPS, bCoup));
			double dn = enuIPS[1] - enuCoup[1];
			double de = enuIPS[0] - enuCoup[0];
			double du = enuIPS[2] - enuCoup[2];
			//fprintf(fpnorm, "%8.1f,%6.3f,%d\n", IPSdata[0],sqrt(dn*dn+de*de+du*du),status);
			fprintf(fpnorm, "%8.1f,%6.3f,%d\n", IPSdata[0], Coupdata[3]-12.6, status);
		}
	}
}


bool ExtractGNGGA(const std::string& inputFile, const std::string& outputFile) {
	std::ifstream in(inputFile.c_str());
	if (!in.is_open()) return false;

	std::ofstream out(outputFile.c_str());
	if (!out.is_open()) {
		in.close();
		return false;
	}

	const std::string target = "";//"THZK RTK_NMEA";
	std::string line;

	while (std::getline(in, line)) {
		size_t pos = line.find(target);
		if (pos != std::string::npos) {
			size_t dataStart = pos + target.length();

			// 跳过空格/制表符
			while (dataStart < line.size() &&
				(line[dataStart] == ' ' || line[dataStart] == '\t')) {
				++dataStart;
			}

			// 验证数据头
			if (line.compare(dataStart, 6, "$GNGGA") == 0) {
				out << line.substr(dataStart) << '\n';
			}
		}
	}

	in.close();
	out.close();
	return true;
}


void extraDaoyYuan(const std::string &filePath)
{
	std::ifstream iFile(filePath);
	if (!iFile.is_open()) {
		throw std::runtime_error("Failed to open input file");
	}
	UTCTIME utct; GPSTIME gpst;
	Matrix Pos;
	double Data[9]{}, XYZ[3]{}, BLH[3];
	std::string line;
	double lever1[] = { 0.385,1.317,-1.338 };    // INS-天线 1.317,0.385,-1.338
	double lever2[] = { 0.028 ,-0.077 ,-0.052 };  //INS-后轮中心 -0.077 ,0.028 ,0.052

	while (1)
	{
		if (!std::getline(iFile, line))break;
		static int skipline = 1;
		if (skipline > 0) {
			skipline--;
			continue;
		}
		std::vector<std::string> tokens = FILEOpt::splitString(line, ',');
		
		//Data[0] = std::stod(tokens[1]);
		//Data[1] = std::stod(tokens[2]); //* Rad;  //纬度
		//Data[2] = std::stod(tokens[3]);// *Rad;  //经度
		//Data[3] = std::stod(tokens[4]);

		//Data[4] = std::stod(tokens[25]);  //航向，俯仰，横滚   输入：横滚，俯仰，航向
		//Data[5] = std::stod(tokens[25+1]);  // 俯仰
		//Data[6] = std::stod(tokens[25+2]);  // 横滚
		//ins_to_gnss(Data + 1, Data + 4, BLH);

		Data[0] = std::stod(tokens[0]);
		Data[1] = std::stod(tokens[2]); //* Rad;  //纬度
		Data[2] = std::stod(tokens[1]);// *Rad;  //经度
		Data[3] = std::stod(tokens[3]);

		Data[4] = std::stod(tokens[7]);  //航向，俯仰，横滚   输入：横滚，俯仰，航向
		Data[5] = std::stod(tokens[8]);  // 俯仰
		Data[6] = std::stod(tokens[9]);  // 横滚
		gnss_to_ins(Data + 1, Data + 4,lever2,BLH);
		ins_to_gnss(BLH, Data + 4,lever1 ,BLH);
		//gnss_to_ins(BLH, Data + 4, BLH);

		static FILE* fptem = NULL;
		if (fptem == NULL) {
			fopen_s(&fptem, "F:\\CDTH\\2025@Couple\\20250310@Desaysv@rtk@dr@Tunnel\\res\\LB_GNSS.txt", "w");
		}
		else
		{
			fprintf(fptem, "%10.3f,%12.7f,%11.7f,%6.3f\n", Data[0], BLH[0], BLH[1], BLH[2]);
		}

		//Pos=INS2GNSS(XYZ, Data + 5);
		//XYZ[0] = Pos(0, 0);
		//XYZ[1] = Pos(1, 0);
		//XYZ[2] = Pos(2, 0);
		//XYZ2LLH(XYZ, BLH);
		printf("%13.9f,%12.9f,%7.4f\n", BLH[1], BLH[0], BLH[2]);
		

	}
}


// 返回固定率
static double getNEMAFixedRatio(std::string filePath,int &len) {
	std::ifstream iFile(filePath);
	if (!iFile.is_open()) {
		throw std::runtime_error("Failed to open input file");
	}
	std::string line, token;
	unsigned long length = 0;
	unsigned long Fixcount = 0;
	while (1) {
		if (!std::getline(iFile, line))break;
		// Check if the line starts with "$GNGGA"
		if (line.find("$GNGGA") != 0)continue;
		std::stringstream ss(line);
		std::vector<std::string> tokens;
		while (std::getline(ss, token, ',')) {
			tokens.push_back(token);
		}
		length++;
		if (std::stoi(tokens[6]) == 4)Fixcount++;
	}
	len = length;
	if (Fixcount == 0)return 0.0;
	return static_cast<double>(Fixcount) / length;
}


void ComparePOS(){
	std::string rootfolderPath= "F:\\CDTH\\202406@CDTH@CropperGNSS\\Lying\\20250102@TH1100@BSDVRS@UGV@ThreeWalls";
	FILEOpt::Filesuffix = ".pos";
	std::vector<std::string> CfileA;
	std::vector<std::string> CfileB;
	std::ifstream ifileA;
	std::ifstream ifileB;

	std::string Compareidx_A = "V5_POSRTKDLLV34TmpV1SNR@EL@Var@COM@POSRTKPOS";  
    std::string Compareidx_B = "V4_POSRTKDLLV34TmpV1SNR@EL@Var@COM@POSRTKPOS"; //POSRTKDLLV337TmpPOS@SNR@EL@Var@Resi

	std::vector<std::string>folderPaths = {
		//"20240912@TH1100_TH1200_UM960@QXVRS@UGV@ThreeWalls",
		//"20240913@TH1100_TH1200_UM960@QXVRS@UGV@TwoWalls",
		//"20241016@TH1100@BSDVRS@UGV@TwoWalls",
		//"20241018@TH1100@BSDVRS@UGV@TwoWalls",
		//"20241128@TH1100_X5@Car@Positec@Opensky",
		//"20241204@TH1100@BSDVRS@UGV@ThreeWalls",
		//"20241203@TH1100@BSDVRS@UGV@Trees",
		//"20250305@TH1100@BSDVRS@UGV@OneWall",
		"最优次优差距优化",
	};
	double Fixed_Change_limit = 1.0;  //规定固定率只有大于该值，才认为有明显变化
	static FILE* fpCompare = nullptr;
	string CompareFilename = Compareidx_A + "_" + Compareidx_B+".txt";
	if (fpCompare==nullptr) {
		fopen_s(&fpCompare, CompareFilename.c_str(), "w");
		fprintf(fpCompare, "%s,%s,%s\n", "$GNCOMPARE", Compareidx_A.c_str(), Compareidx_B.c_str());
	}
	double countALLFile[2] = {0};
	int All_A2B_change[3] = { 0 };
	for (auto i = folderPaths.begin(); i !=folderPaths.end(); i++)
	{
		printf("%s", "正在处理文件夹 ");
		std::cout << *i << std::endl;
		CfileA = FILEOpt::traversePdir(rootfolderPath + "\\" + *i + "\\" + Compareidx_A);
		CfileB = FILEOpt::traversePdir(rootfolderPath + "\\" + *i + "\\" + Compareidx_B);

		vector<double>FixedCount[2];
		int A2B_change[3] = { 0 }, lenth[2]{};//变化量和历元数，历元数未启用
		for (size_t j = 0; j < CfileA.size(); j++)
		{
			FixedCount[0].push_back(getNEMAFixedRatio(CfileA[j],lenth[0]) * 100.0);
			FixedCount[1].push_back(getNEMAFixedRatio(CfileB[j],lenth[1]) * 100.0);
			if ((FixedCount[0][j] - FixedCount[1][j])> Fixed_Change_limit)A2B_change[0]++;   // A比B版本有多少个文件有改善
			else if (fabs(FixedCount[0][j] - FixedCount[1][j]) < Fixed_Change_limit)A2B_change[1]++;   // A比B有多少个文件无变化
			else  A2B_change[2]++;   // B比A有多少个文件固定率下降
			if (fpCompare)
			{
				fprintf(fpCompare, "%s,%5.2f,%5.2f,%6.2f,%d\n", (*i).c_str(), FixedCount[0][j], FixedCount[1][j], FixedCount[0][j]- FixedCount[1][j],j+1);
				printf("%s,%5.2f,%5.2f,%6.2f,%d\n", (*i).c_str(), FixedCount[0][j], FixedCount[1][j], FixedCount[0][j] - FixedCount[1][j],j+1);
			}
		}
		fprintf(fpCompare, "%s,%5.2f,%5.2f\n", "$SUM 该文件夹平均固定率", mean(FixedCount[0]), mean(FixedCount[1]));
		fprintf(fpCompare, "%s%d,%d,%d,%s\n","$SUM 相对于B版本，有：",A2B_change[0],A2B_change[1],A2B_change[2],"上升、不变、下降");
		printf("%s,%d,%d,%d,%s\n", "相对于B版本，有：", A2B_change[0], A2B_change[1], A2B_change[2], "上升、不变、下降");
		//fprintf(fpCompare,"%s,%5.2f,%5.2f, 相对于B文件有 %d%s, 总共有 %d 个文件\n", "$SUM 该文件夹平均固定率", mean(FixedCount[0]), mean(FixedCount[1]), A2B_change, "个文件固定率上升", CfileA.size());
		//printf("%s,%5.2f,%5.2f, 相对于B文件有 %d%s, 总共有 %d 个文件\n", "$SUM 该文件夹平均固定率", mean(FixedCount[0]), mean(FixedCount[1]), A2B_change, "个文件固定率上升", CfileA.size());
	
		countALLFile[0] += mean(FixedCount[0]);
		countALLFile[1] += mean(FixedCount[1]);
		for (int k = 0; k < 3;k++)All_A2B_change[k] += A2B_change[k];
	}
	countALLFile[0] = countALLFile[0] / folderPaths.size();
	countALLFile[1] = countALLFile[1] / folderPaths.size();

	fprintf(fpCompare, "\n%s,%5.2f,%5.2f\n", "$ALL 所有文件夹中平均固定率", countALLFile[0], countALLFile[1]);
	fprintf(fpCompare, "%s%d,%d,%d,%s\n", "$ALL 相对于B版本，有：", All_A2B_change[0], All_A2B_change[1], All_A2B_change[2], "上升、不变、下降");
	printf("\n%s,%5.2f,%5.2f\n", "$ALL 所有文件夹中平均固定率", countALLFile[0], countALLFile[1]);
	printf("%s%d,%d,%d,%s\n", "$ALL 相对于B版本，有：", All_A2B_change[0], All_A2B_change[1], All_A2B_change[2], "上升、不变、下降");
	
	std::cout << "遍历文件夹数量   " << FILEOpt::Foldercount << std::endl;
	std::cout << "遍历文件数量    " << FILEOpt::Filecount << std::endl;
}





bool DecodeNEMA(int PosType)
{
	std::string folderPath;
	std::cout << "请输入文件路径:定位的pos文件夹路径)" << std::endl << std::endl;
	std::getline(std::cin, folderPath);  // 从控制台读取文件夹路径

	std::ifstream inputFile(folderPath);
	if (!inputFile.is_open()) {
		throw std::runtime_error("Failed to open input file: " + folderPath);
	}

	std::string outputPath = FILEOpt::changeStrExtension(folderPath, ".txt");
	UTCTIME utc; GPSTIME gpsNEMA;
	std::string line;
	std::string token;
	double XYZ[3]{}, BLH[3]{};
	int solStatus = 0;
	std::vector<double>Pos[3];
	double Rad = PI / 180;
	double DaoYuanPos[9] = { 0.0 };
	while (std::getline(inputFile, line)) {
		// Check if the line starts with "$GNGGA"    $ACCUR
		if (line.find("$ACCUR") == 0 && PosType == 0) {
			std::vector<std::string>tokens = FILEOpt::splitString(line, ',');
			utc.Year = std::stoi(tokens[4].substr(0, 4));
			utc.Month = std::stoi(tokens[4].substr(4, 2));
			utc.Day = std::stoi(tokens[4].substr(6, 2));
		}
		if (line.find("$GNGGA") == 0 && PosType == 0) {

			std::stringstream ss(line);
			std::vector<std::string> tokens;
			while (std::getline(ss, token, ',')) {
				tokens.push_back(token);
			}
			if (1) {
				BLH[0] = std::stod(tokens[2]) / 100.0;
				BLH[1] = std::stod(tokens[4]) / 100.0;
				BLH[2] = std::stod(tokens[9]);
				for (size_t i = 0; i < 2; i++)BLH[i] = DegMinSec2D(BLH[i]);
				utc.Hour = std::stoi(tokens[1].substr(0, 2));
				utc.Minute = std::stoi(tokens[1].substr(2, 2));
				utc.Second = std::stod(tokens[1].substr(4));
				solStatus = std::stoi(tokens[6]);
				gpsNEMA = UTCtoGPST(utc);

				static FILE* Posfp = NULL;
				if (!Posfp)fopen_s(&Posfp, outputPath.c_str(), "w");
				if (gpsNEMA.fracSec < 0.01)fprintf(Posfp, "%4d %8.2f %10.6f %11.6f %6.2f %d\n", gpsNEMA.week, gpsNEMA.secWeek + gpsNEMA.fracSec, BLH[0], BLH[1], BLH[2], solStatus);
			}
			if (tokens[6] == "4")
			{
				BLH[0] = std::stod(tokens[2]) / 100.0;
				BLH[1] = std::stod(tokens[4]) / 100.0;
				BLH[2] = std::stod(tokens[9]);
				for (size_t i = 0; i < 2; i++)BLH[i] = DegMinSec2D(BLH[i]) * Rad;  //度分转为度

				if (tokens[3] == "S")BLH[0] = -BLH[0];
				if (tokens[5] == "W")BLH[1] = -BLH[1];
				LLH2XYZ(BLH, XYZ);
				for (size_t i = 0; i < 3; i++)Pos[i].push_back(XYZ[i]);
			}

		}
		if (PosType == 1) {
			
			static FILE* Posfp = NULL;
			if (!Posfp) {
				fopen_s(&Posfp, outputPath.c_str(), "w");
				std::getline(inputFile, line);
			}
			extraDaoYuanPos(line, DaoYuanPos);
			if (1) {
				
				fprintf(Posfp, "%8.2f %10.6f %11.6f %6.2f\n", DaoYuanPos[0], DaoYuanPos[1], DaoYuanPos[2], DaoYuanPos[3]);
			}
			for (int i = 0; i < 9; i++)
			{
				std::getline(inputFile, line); // 采样率100Hz，转为1Hz输出
			}
			//// IE 的HW格式
			//if (line.find("Fixed") != std::string::npos) {
			//	std::stringstream ss(line);
			//	std::vector<std::string> tokens;
			//	tokens = FILEOpt::splitString(line, ' ');
			//	if (tokens[5] != "1")continue;
			//	for (size_t i = 0; i < 3; i++)XYZ[i] = stod(tokens[16 + i]);
			//	for (size_t i = 0; i < 3; i++)Pos[i].push_back(XYZ[i]);
			//}
		}

		if (PosType == 2) {
			if (line.find("#INSPVAXA") == std::string::npos)continue;
			static FILE* Posfp = NULL;
			if (!Posfp) {
				fopen_s(&Posfp, outputPath.c_str(), "w");
			}
			extraDeSayRealPos(line, DaoYuanPos);
			fprintf(Posfp, "%8.2f %10.6f %11.6f %6.2f\n", DaoYuanPos[0], DaoYuanPos[1], DaoYuanPos[2], DaoYuanPos[3]);
			for (int i = 0; i < 99; i++)
			{
				std::getline(inputFile, line); // 采样率100Hz，转为1Hz输出
			}
		}
	}
	if (Pos[0].size() == 0)return false;

	double RefPos[3]{};
	for (int j = 0; j < 3; j++)RefPos[j] = meanWithoutOutliers(Pos[j]);
	//for (size_t i = 0; i < Pos[0].size(); i++)
	//{
	//    for (int j = 0; j < 3; j++)RefPos[j] += Pos[j][i];
	//}
	//for (int j = 0; j < 3; j++)RefPos[j] = RefPos[j]/ Pos[0].size();

	std::cout << "该文件固定解的平均值为 XYZ： " << std::endl;
	printf("%12.4f,%12.4f,%12.4f\n", RefPos[0], RefPos[1], RefPos[2]);
	XYZ2LLH(RefPos, BLH);

	std::cout << "该文件固定解的平均值为 BLH： " << std::endl;
	printf("%12.9f,%13.9f,%8.4f\n", BLH[0] * 180 / PI, BLH[1] * 180 / PI, BLH[2]);

	printf("\n\n");
	return true;
}




//Test
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#define NUM_TESTS 1000

constexpr const int MAXOBSLIM = 50;
void generate_positive_definite_matrix(int n, double* matrix) {

	double L[MAXOBSLIM * MAXOBSLIM] = { 0.0 };

	// 生成下三角矩阵 L，对角线元素为随机正数，非对角线元素为随机数
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i < j) {
				L[i * n + j] = 0.0;  // 上三角部分为零
			}
			else {
				double value = (double)rand() / RAND_MAX;  // [0, 1] 范围内的随机数
				if (i == j) {
					L[i * n + j] = value + 0.2;  // 对角线元素增加0.1以确保正定性
				}
				else {
					L[i * n + j] = value * 2.0 - 1.0;  // 下三角非对角线元素为 [-1, 1]
				}
			}
		}
	}

	// 计算 M = L * L^T，确保正定对称
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < n; k++) {
				sum += L[i * n + k] * L[j * n + k];  // 计算点积
			}
			matrix[i * n + j] = sum;
		}
	}
}

#define MAXLIM MAXOBSLIM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXLIM 100
#define EPS 1e-10

// 交换矩阵中的两行
void swap_rows(double mat[], int i, int j, int n) {
	for (int k = 0; k < n; k++) {
		double temp = mat[i * n + k];
		mat[i * n + k] = mat[j * n + k];
		mat[j * n + k] = temp;
	}
}

// 计算整数 gcd
static int gcd_int(int a, int b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	while (b != 0) {
		int t = a % b;
		a = b;
		b = t;
	}
	return a;
}

// 对矩阵进行初等行变换
static bool transform_matrix(int n, const double Z[], double B[], double P[]) {
	// 复制原始矩阵Z到B
	for (int i = 0; i < n * n; i++) {
		B[i] = Z[i];
	}

	// 初始化变换矩阵P为单位矩阵
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			P[i * n + j] = (i == j) ? 1.0 : 0.0;
		}
	}

	// 最后一列的索引
	const int last_col = n - 1;
	const int last_row = n - 1;

	// -------- gcd 检查 --------
	int g = 0;
	for (int i = 0; i < n; i++) {
		int val = (int)round(B[i * n + last_col]);
		g = gcd_int(g, val);
	}
	if (g != 1) {
		// gcd 不是 1，无法规约成 ±1，直接返回
		return false;
	}

	// 步骤1：确保最后一行最后一列元素为±1
	bool found_unit = false;
	while (!found_unit) {
		// 检查最后一行最后一列是否为±1
		if (fabs(fabs(B[last_row * n + last_col]) - 1.0) < EPS) {
			found_unit = true;
			continue;
		}

		// 检查其他行是否有±1
		int pivot_row = -1;
		for (int i = 0; i < n - 1; i++) {
			if (fabs(fabs(B[i * n + last_col]) - 1.0) < EPS) {
				pivot_row = i;
				break;
			}
		}

		if (pivot_row != -1) {
			// 交换行到最底部
			swap_rows(B, pivot_row, last_row, n);
			swap_rows(P, pivot_row, last_row, n);
			found_unit = true;
		}
		else {
			// 寻找最小非零元素
			int min_row = -1;
			double min_val = 0.0;
			for (int i = 0; i < n; i++) {
				double abs_val = fabs(B[i * n + last_col]);
				if (abs_val > EPS) {
					if (min_row == -1 || abs_val < min_val) {
						min_row = i;
						min_val = abs_val;
					}
				}
			}

			if (min_row == -1) {
				// 最后一列全为零（不可能发生，因为矩阵可逆）
				break;
			}

			// 将最小元素所在行交换到最后
			if (min_row != last_row) {
				swap_rows(B, min_row, last_row, n);
				swap_rows(P, min_row, last_row, n);
			}

			// 使用最后一行消去其他行
			for (int i = 0; i < n - 1; i++) {
				if (fabs(B[i * n + last_col]) > EPS) {
					double factor = round(B[i * n + last_col] / B[last_row * n + last_col]);
					for (int j = 0; j < n; j++) {
						B[i * n + j] -= factor * B[last_row * n + j];
						P[i * n + j] -= factor * P[last_row * n + j];
					}
				}
			}
		}
	}

	// 步骤2：确保最后一行最后一列元素为+1
	if (fabs(B[last_row * n + last_col] + 1.0) < EPS) {
		// 如果是-1，则乘以-1
		for (int j = 0; j < n; j++) {
			B[last_row * n + j] *= -1;
			P[last_row * n + j] *= -1;
		}
	}

	// 步骤3：使用最后一行消去其他行的最后一列
	for (int i = 0; i < n - 1; i++) {
		if (fabs(B[i * n + last_col]) > EPS) {
			double factor = round(B[i * n + last_col] / B[last_row * n + last_col]);
			for (int j = 0; j < n; j++) {
				B[i * n + j] -= factor * B[last_row * n + j];
				P[i * n + j] -= factor * P[last_row * n + j];
			}
		}
	}

	// 清除浮点误差（将接近零的值设为零）
	for (int i = 0; i < n - 1; i++) {
		if (fabs(B[i * n + last_col]) < EPS) {
			B[i * n + last_col] = 0.0;
		}
	}
	return true;
}


// 计算矩阵行列式的函数
// 参数：matrix - 一维数组表示的矩阵（按行优先存储）
//        n     - 矩阵的阶数
// 返回值：矩阵的行列式值
double determinant(double matrix[], int n) {
	// 检查矩阵是否为空或超过最大空间
	if (n <= 0 || n * n > MAXOBSLIM*MAXOBSLIM) {
		printf("Invalid matrix size.\n");
		return 0.0;
	}

	// 创建临时矩阵进行运算（避免修改原矩阵）
	double temp[MAXOBSLIM* MAXOBSLIM];
	for (int i = 0; i < n * n; i++) {
		temp[i] = matrix[i];
	}

	double det = 1.0; // 初始化行列式值
	int sign = 1;     // 用于记录行交换的次数

	// 高斯消元过程
	for (int col = 0; col < n; col++) {
		// 寻找主元（当前列中绝对值最大的元素）
		int max_row = col;
		double max_val = fabs(temp[col * n + col]);

		for (int row = col + 1; row < n; row++) {
			double element = fabs(temp[row * n + col]);
			if (element > max_val) {
				max_val = element;
				max_row = row;
			}
		}

		// 如果主元为0，则行列式为0
		if (max_val < 1e-9) {
			return 0.0;
		}

		// 如果需要，交换行
		if (max_row != col) {
			sign *= -1; // 行交换改变行列式符号
			for (int j = 0; j < n; j++) {
				double swap = temp[col * n + j];
				temp[col * n + j] = temp[max_row * n + j];
				temp[max_row * n + j] = swap;
			}
		}

		// 消元：使当前列下方的元素变为0
		for (int row = col + 1; row < n; row++) {
			double factor = temp[row * n + col] / temp[col * n + col];
			for (int j = col; j < n; j++) {
				temp[row * n + j] -= factor * temp[col * n + j];
			}
		}

		// 累乘主对角线元素
		det *= temp[col * n + col];
	}

	return det * sign; // 应用行交换的符号
}

void test_incremental_vs_direct() {
	srand(42);  // 确保结果可重复

	// 生成随机正定矩阵和浮点解
	double Qahat[MAXOBSLIM * MAXOBSLIM];
	double ahat[MAXOBSLIM];


	// 为两种方法准备工作空间
	double Qzhat_inc[MAXOBSLIM * MAXOBSLIM], Z_inc[MAXOBSLIM * MAXOBSLIM];
	double L_inc[MAXOBSLIM * MAXOBSLIM], D_inc[MAXOBSLIM];
	double zhat_inc[MAXOBSLIM], iZt_inc[MAXOBSLIM * MAXOBSLIM];

	double Qzhat_dir[MAXOBSLIM * MAXOBSLIM], Z_dir[MAXOBSLIM * MAXOBSLIM];
	double L_dir[MAXOBSLIM * MAXOBSLIM], D_dir[MAXOBSLIM];
	double zhat_dir[MAXOBSLIM], iZt_dir[MAXOBSLIM * MAXOBSLIM];

	// 创建降维后的矩阵和向量
	double Qrr[(MAXOBSLIM - 1) * (MAXOBSLIM - 1)];
	double ahat_rr[MAXOBSLIM - 1];
	double Qtt[(MAXOBSLIM - 1) * (MAXOBSLIM - 1)] = {0.0};
	double ZT[MAXOBSLIM * MAXOBSLIM] = { 0.0 };

	double start, end;

	// 测试维度范围
	int dims[] = { 15, 20, 25, 30,35,40 };   //别超过MAXOBSLIM
	int num_dims = sizeof(dims) / sizeof(dims[0]);
	int nIncremFail[6] = {0};
	double B[MAXOBSLIM * MAXOBSLIM], P[MAXOBSLIM * MAXOBSLIM];
	double B11[MAXOBSLIM * MAXOBSLIM] = { 0.0 };
	bool isSuc = false;

	// 结果存储
	double incremental_times[MAXOBSLIM];
	double direct_times[MAXOBSLIM];

	// 设置参数并复制子矩阵
	int s1[2] = { 0, 0 };
	int s2[2] = { 0, 0 };
	double Qbb[MAXOBSLIM * MAXOBSLIM];

	

	for (int d = 0; d < num_dims; d++) {
		int n = dims[d];
		printf("Testing dimension: %d\n", n);

		clock_t incremental_total = 0;
		clock_t direct_total = 0;

		for (int test = 0; test < NUM_TESTS; test++) {

			generate_positive_definite_matrix(n, Qahat);
			for (int i = 0; i < n; i++) {
				ahat[i] = (double)rand() / RAND_MAX;
			}

			// 首先进行完整的降相关作为起点
			if (decorrel(Qahat, ahat, Qzhat_inc, Z_inc, L_inc, D_inc, zhat_inc, iZt_inc, n) == false) {
				nIncremFail[d]++;
				continue;  // 不比较
			}
			MatrixTranspose1(n, n, Z_inc, ZT);   // 维度别超过MAXOBSLIM
			MatrixCopySub(s1,s2, n - 1, n - 1, n, n, Qahat, n - 1, n - 1, Qrr, false);
			if (1) {
				// 测试提供初值的增量更新方法
				start = clock();
				isSuc = transform_matrix(n, ZT, B, P);   // B=PZ，其中P为初等行变换矩阵，B的前n-1行n-1列也是整数单模矩阵
				if (!isSuc) {
					nIncremFail[d]++;
					continue;  // 不比较
				}
				MatrixMultiply_HPHT(n, n, B, Qahat, false, Qtt);
				MatrixCopySub(s1, s2, n - 1, n - 1, n, n, Qtt, n - 1, n - 1, Qbb, false);
				decorrel(Qbb, ahat, Qzhat_dir, Z_dir, L_dir, D_dir, zhat_dir, iZt_dir, n - 1);
				MatrixTranspose1(n - 1, n - 1, Z_dir, ZT);   // 维度别超过MAXOBSLIM
				MatrixCopySub(s1, s2, n - 1, n - 1, n, n, B, n - 1, n - 1, B11, false);
				MatrixMultiply(n - 1, n - 1, ZT, n - 1, n - 1, B11, Z_dir);  //通过增量形式获得新的Z转换矩阵

				if (0) {
					MatrixMultiply_HPHT(n - 1, n - 1, Z_dir, Qrr, false, Qzhat_dir);
					PrintMatrix(Qzhat_dir, n - 1, n - 1);
				}


				end = clock();
				incremental_total += (end - start);
			}

			if (1) {
				// 测试对Qrr从头分解的方法
				start = clock();
				decorrel(Qrr, ahat, Qzhat_dir, Z_dir, L_dir, D_dir, zhat_dir, iZt_dir, n - 1);
				//PrintMatrix(Qzhat_dir, n - 1, n - 1);
				end = clock();
				direct_total += (end - start);
			}


			
		}

		// 计算平均时间（毫秒）
		incremental_times[d] = ((double)incremental_total / CLOCKS_PER_SEC) * 1000 / NUM_TESTS;
		direct_times[d] = ((double)direct_total / CLOCKS_PER_SEC) * 1000 / NUM_TESTS;

		printf("  Incremental avg: %.4f ms\n", incremental_times[d]);
		printf("  Direct avg: %.4f ms\n", direct_times[d]);
		printf("  Speedup: %.2fx\n", direct_times[d] / incremental_times[d]);
	}

	// 打印汇总结果
	printf("\n===== Summary =====\n");
	printf("Dim\tIncremental(ms)\tDirect(ms)\tSpeedup\n");
	for (int d = 0; d < num_dims; d++) {
		printf("%d\t%.4f\t\t%.4f\t\t%.2fx\n",
			dims[d], incremental_times[d], direct_times[d],
			direct_times[d] / incremental_times[d]);

		//printf("%d%s%d\n", dims[d], "维增量降相关失败次数:", nIncremFail[d]);
	}

}


void test_multi_removal() {
	srand(42);  // 确保结果可重复

	// 生成随机正定矩阵和浮点解
	double Qahat[MAXOBSLIM * MAXOBSLIM];
	double ahat[MAXOBSLIM];


	// 为两种方法准备工作空间
	double Qzhat_inc[MAXOBSLIM * MAXOBSLIM], Z_inc[MAXOBSLIM * MAXOBSLIM];
	double L_inc[MAXOBSLIM * MAXOBSLIM], D_inc[MAXOBSLIM];
	double zhat_inc[MAXOBSLIM], iZt_inc[MAXOBSLIM * MAXOBSLIM];

	double Qzhat_dir[MAXOBSLIM * MAXOBSLIM], Z_dir[MAXOBSLIM * MAXOBSLIM];
	double L_dir[MAXOBSLIM * MAXOBSLIM], D_dir[MAXOBSLIM];
	double zhat_dir[MAXOBSLIM], iZt_dir[MAXOBSLIM * MAXOBSLIM];

	// 创建降维后的矩阵和向量
	double Qrr[(MAXOBSLIM - 1) * (MAXOBSLIM - 1)];
	double ahat_rr[MAXOBSLIM - 1];
	double Qtt[(MAXOBSLIM - 1) * (MAXOBSLIM - 1)] = { 0.0 };
	double ZT[MAXOBSLIM * MAXOBSLIM] = { 0.0 };

	double start, end;

	// 测试维度范围
	int dims[] = { 15, 20, 25, 30,35,40 };   //别超过MAXOBSLIM
	int num_dims = sizeof(dims) / sizeof(dims[0]);
	int ndel = 10;
	int nIncremFail[6] = { 0 };


	double B[MAXOBSLIM * MAXOBSLIM], P[MAXOBSLIM * MAXOBSLIM];
	double B11[MAXOBSLIM * MAXOBSLIM] = { 0.0 };
	bool isSuc = false;

	// 结果存储
	double incremental_times[MAXOBSLIM];
	double direct_times[MAXOBSLIM];

	// 设置参数并复制子矩阵
	int s1[2] = { 0, 0 };
	int s2[2] = { 0, 0 };
	double Qbb[MAXOBSLIM * MAXOBSLIM];



	for (int d = 0; d < num_dims; d++) {
		int n = dims[d];
		printf("Testing dimension: %d\n", n);

		clock_t incremental_total = 0;
		clock_t direct_total = 0;

		for (int test = 0; test < NUM_TESTS; test++) {

			generate_positive_definite_matrix(n, Qahat);
			for (int i = 0; i < n; i++) {
				ahat[i] = (double)rand() / RAND_MAX;
			}


			for (size_t del = 0; del < ndel; del++)
			{
				int row = n - del;

				
				// 首先进行完整的降相关作为起点
				if (decorrel(Qahat, ahat, Qzhat_inc, Z_inc, L_inc, D_inc, zhat_inc, iZt_inc, row) == false) {
					nIncremFail[d]++;
					break;  // 不比较
				}
				MatrixTranspose1(row, row, Z_inc, ZT);   // 维度别超过MAXOBSLIM
				MatrixCopySub(s1, s2, row - 1, row - 1, row, row, Qahat, row - 1, row - 1, Qrr, false);
				if (1) {
					// 测试提供初值的增量更新方法
					start = clock();
					isSuc = transform_matrix(row, ZT, B, P);   // B=PZ，其中P为初等行变换矩阵，B的前row-1行row-1列也是整数单模矩阵
					if (!isSuc) {
						nIncremFail[d]++;
						continue;  // 不比较
					}
					MatrixMultiply_HPHT(row, row, B, Qahat, false, Qtt);
					MatrixCopySub(s1, s2, row - 1, row - 1, row, row, Qtt, row - 1, row - 1, Qbb, false);
					decorrel(Qbb, ahat, Qzhat_dir, Z_dir, L_dir, D_dir, zhat_dir, iZt_dir, row - 1);
					MatrixTranspose1(row - 1, row - 1, Z_dir, ZT);   // 维度别超过MAXOBSLIM
					MatrixCopySub(s1, s2, row - 1, row - 1, row, row, B, row - 1, row - 1, B11, false);
					MatrixMultiply(row - 1, row - 1, ZT, row - 1, row - 1, B11, Z_dir);  //通过增量形式获得新的Z转换矩阵

					if (0) {
						MatrixMultiply_HPHT(row - 1, row - 1, Z_dir, Qrr, false, Qzhat_dir);
						PrintMatrix(Qzhat_dir, row - 1, row - 1);
					}


					end = clock();
					incremental_total += (end - start);
				}

				if (1) {
					// 测试对Qrr从头分解的方法

					memset(Qzhat_dir, 0.0, sizeof(Qzhat_dir));
					memset(Z_dir, 0.0, sizeof(Z_dir));
					memset(L_dir, 0.0, sizeof(L_dir));
					memset(D_dir, 0.0, sizeof(D_dir));
					memset(iZt_dir, 0.0, sizeof(iZt_dir));


					start = clock();
					decorrel(Qrr, ahat, Qzhat_dir, Z_dir, L_dir, D_dir, zhat_dir, iZt_dir, row - 1);
					//PrintMatrix(Qzhat_dir, row - 1, row - 1);
					end = clock();
					direct_total += (end - start);
				}

				memset(Qzhat_dir, 0.0, sizeof(Qzhat_dir));
				memset(Z_dir, 0.0, sizeof(Z_dir));
				memset(L_dir, 0.0, sizeof(L_dir));
				memset(D_dir, 0.0, sizeof(D_dir));
				memset(iZt_dir, 0.0, sizeof(iZt_dir));
				memset(ZT, 0.0, sizeof(ZT));
				memset(B, 0.0, sizeof(B)); 
				memset(Qbb, 0.0, sizeof(Qbb));
				memset(B11, 0.0, sizeof(B11));


				memset(Qzhat_inc, 0.0, sizeof(Qzhat_inc));
				memset(Z_inc, 0.0, sizeof(Z_inc));
				memset(L_inc, 0.0, sizeof(L_inc));
				memset(D_inc, 0.0, sizeof(D_inc));
				memset(iZt_inc, 0.0, sizeof(iZt_inc));
				double Temp[MAXOBSLIM * MAXOBSLIM] = { 0.0 };

				//PrintMatrix(Qahat, row, row);
				MatrixCopySub(s1, s2, row - 1, row - 1, row, row, Qahat, row - 1, row - 1, Temp, false);
				//PrintMatrix(Temp,row-1,row-1);
				for (size_t i = 0; i < MAXOBSLIM*MAXOBSLIM; i++)
				{
					Qahat[i] = Temp[i];
				}
				//PrintMatrix(Qahat, row - 1, row - 1);
			}

			



		}

		// 计算平均时间（毫秒）
		incremental_times[d] = ((double)incremental_total / CLOCKS_PER_SEC) * 1000 / NUM_TESTS;
		direct_times[d] = ((double)direct_total / CLOCKS_PER_SEC) * 1000 / NUM_TESTS;

		printf("  Incremental avg: %.4f ms\n", incremental_times[d]);
		printf("  Direct avg: %.4f ms\n", direct_times[d]);
		printf("  Speedup: %.2fx\n", direct_times[d] / incremental_times[d]);
	}

	// 打印汇总结果
	printf("\n===== Summary =====\n");
	printf("Dim\tIncremental(ms)\tDirect(ms)\tSpeedup\n");
	for (int d = 0; d < num_dims; d++) {
		printf("%d\t%.4f\t\t%.4f\t\t%.2fx\n",
			dims[d], incremental_times[d], direct_times[d],
			direct_times[d] / incremental_times[d]);

		printf("%d%s%d\n", dims[d], "维增量降相关失败次数:", nIncremFail[d]);
	}
}