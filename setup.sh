#!/bin/bash

# For MacOS and Linux.
# Run `source install.sh`

# 检查是否安装了 Conda
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed. Please install Miniconda or Anaconda first."
    exit 1
fi

# 使 `conda activate` 在 shell 中生效
source $(conda info --base)/etc/profile.d/conda.sh

# 初始化 Conda（可选，仅第一次运行需要）
conda init

# 创建并激活 Conda 虚拟环境
ENV_NAME="h9-env"
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
else
    conda create -n "$ENV_NAME" python=3.11 -y
fi

conda activate "$ENV_NAME"

# 卸载冲突的 Pandas 和 NumPy 版本
python -m pip uninstall pandas -y
python -m pip uninstall numpy -y

# 安装 Biopython, Pandas 和 MAFFT
conda install -c bioconda -c conda-forge biopython=1.84 pandas mafft -y

# 获取 Conda 环境的 `bin` 目录
CONDA_BIN_DIR="${CONDA_PREFIX}/bin"

# 获取操作系统和架构信息
OS=$(uname -s)
ARCH=$(uname -m)

# 匹配 BLAST 二进制文件
BASE="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/"
case "$OS-$ARCH" in
    Linux-x86_64) BLAST_URL="ncbi-blast-2.16.0+-x64-linux.tar.gz" ;;
    Linux-aarch64) BLAST_URL="ncbi-blast-2.16.0+-aarch64-linux.tar.gz" ;;
    Darwin-x86_64) BLAST_URL="ncbi-blast-2.16.0+-x64-macosx.tar.gz" ;;
    Darwin-arm64|Darwin-aarch64) BLAST_URL="ncbi-blast-2.16.0+-aarch64-macosx.tar.gz" ;;
    *) echo "Error: Unsupported platform: $OS-$ARCH. Please install BLAST manually."; exit 1 ;;
esac

# 下载并解压 BLAST
echo "Installing BLAST from $BASE$BLAST_URL"
curl -O "${BASE}${BLAST_URL}"

# 获取 BLAST 文件夹名
BLAST_FOLDER=$(tar -tzf "${BLAST_URL}" | head -1 | cut -f1 -d"/")

# 解压 BLAST 并移动到 Conda 环境的 bin 目录
tar -xzf "${BLAST_URL}"
rm "${BLAST_URL}"

echo "Moving BLAST binaries to Conda environment: ${CONDA_BIN_DIR}"
mv "${BLAST_FOLDER}/bin/"* "${CONDA_BIN_DIR}/"
rm -rf "${BLAST_FOLDER}"

# 设置环境变量
export H9_PYTHON_PATH="${CONDA_BIN_DIR}/python"
export H9_MAFFT_PATH="${CONDA_BIN_DIR}/mafft"
export H9_BLAST_PATH="${CONDA_BIN_DIR}/blastn"

# 将环境变量写入 .bashrc / .zshrc 以便长期使用
if [ "$SHELL" = "/bin/zsh" ]; then
    RC_FILE="$HOME/.zshrc"
elif [ "$SHELL" = "/bin/bash" ]; then
    RC_FILE="$HOME/.bashrc"
else
    RC_FILE="$HOME/.profile"
fi

echo "Exporting environment variables to $RC_FILE..."
{
    echo "export H9_PYTHON_PATH=${H9_PYTHON_PATH}"
    echo "export H9_MAFFT_PATH=${H9_MAFFT_PATH}"
    echo "export H9_BLAST_PATH=${H9_BLAST_PATH}"
} >> "$RC_FILE"

echo "Installation complete. Please run 'source $RC_FILE' or restart your shell to apply changes."

# 显示安装的版本信息
echo "Checking installed versions..."
echo "Python version: $(python --version)"
echo "MAFFT version: $(mafft --version)"
echo "BLAST version: $($H9_BLAST_PATH -version | head -n 1)"

# 检查安装的 Python 包版本
echo "Checking installed package versions..."
echo -n "Biopython version: " && python -c "import Bio; print(Bio.__version__)"
echo -n "Pandas version: " && python -c "import pandas as pd; print(pd.__version__)"
echo -n "NumPy version: " && python -c "import numpy as np; print(np.__version__)"

echo "Installation and verification complete."