import os

# 设置默认基因和输出文件
gtfFile = "human.gtf"
gene = "FSTL1"

# 检查是否有传入基因名作为命令行参数
import sys
if len(sys.argv) > 1:
    gene = sys.argv[1]

# 打开输出文件
with open("singleGeneExpFSTL1.txt", "w") as wf:
    wf.write("Id\t{}\tType\tCancerType\n".format(gene))

    # 遍历当前目录中的文件
    for file in os.listdir("."):
        if file.startswith("symbol.") and file.endswith(".txt"):
            # 提取 "type" 从文件名
            type = file.split(".")[1]

            with open(file, "r") as rf:
                sampleArr = []
                for line in rf:
                    arr = line.strip().split("\t")
                    if not sampleArr:
                        sampleArr = arr
                    elif gene == arr[0]:
                        for i in range(1, len(arr)):
                            sampleName = sampleArr[i]
                            barcodeArr = sampleName.split("-")
                            if barcodeArr[3].startswith("0"):
                                sampleType = "Tumor"
                            else:
                                sampleType = "Normal"
                            wf.write(f"{sampleName}\t{arr[i]}\t{sampleType}\t{type}\n")
                        break
