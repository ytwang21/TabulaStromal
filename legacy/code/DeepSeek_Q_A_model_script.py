#!/usr/bin/env python3
import os
import glob
import xml.etree.ElementTree as ET
import logging
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
from tqdm import tqdm
import pandas as pd
import re

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DeepSeekProcessor:
    def __init__(self, model_path, max_length=4096):
        """
        初始化DeepSeek模型
        Args:
            model_path: 模型路径
            max_length: 最大生成长度
        """
        self.max_length = max_length
        logger.info(f"正在加载模型: {model_path}")
        self.tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
        self.model = AutoModelForCausalLM.from_pretrained(
            model_path,
            torch_dtype=torch.float16,
            device_map="auto",
            trust_remote_code=True
        )
        # 设置pad_token
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token

        logger.info("模型加载完成")

    def extract_text_from_xml(self, xml_file):
        """
        从XML文件中提取文本内容
        """
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            # 提取所有文本内容
            text_parts = []
            # 递归提取所有元素的文本
            def extract_text(element):
                if element.text and element.text.strip():
                    text_parts.append(element.text.strip())
                for child in element:
                    extract_text(child)
                if element.tail and element.tail.strip():
                    text_parts.append(element.tail.strip())

            extract_text(root)
            # 合并文本并清理
            full_text = ' '.join(text_parts)
            full_text = re.sub(r'\s+', ' ', full_text)  # 合并多余空格

            # 如果文本太长，进行截断
            if len(full_text) > 8000:
                full_text = full_text[:8000] + "... [truncated]"

            return full_text

        except Exception as e:
            logger.error(f"解析XML文件 {xml_file} 时出错: {e}")
            return None

    def ask_question(self, context, question):
        """
        向模型提问并获取回答
        """
        try:
            # 构建提示
            prompt = f"""请基于以下生物医学实验描述文本回答问题。

文本内容:
{context}

问题: {question}

请直接给出答案，不要添加额外的解释。"""

            # 编码输入
            inputs = self.tokenizer(prompt, return_tensors="pt", truncation=True, max_length=4096)

            # 将输入移动到模型所在的设备
            inputs = {k: v.to(self.model.device) for k, v in inputs.items()}

            # 生成回答
            with torch.no_grad():
                outputs = self.model.generate(
                    **inputs,
                    max_new_tokens=512,
                    do_sample=True,
                    temperature=0.7,
                    top_p=0.9,
                    pad_token_id=self.tokenizer.eos_token_id
                )

            # 解码回答
            response = self.tokenizer.decode(outputs[0], skip_special_tokens=True)

            # 提取生成的回答部分（去掉原始提示）
            answer = response[len(prompt):].strip()

            return answer

        except Exception as e:
            logger.error(f"模型推理时出错: {e}")
            return f"Error: {str(e)}"

    def process_single_file(self, xml_file, questions):
        """
        处理单个XML文件
        """
        logger.info(f"处理文件: {os.path.basename(xml_file)}")

        # 提取XML文本内容
        context = self.extract_text_from_xml(xml_file)
        if context is None:
            return None

        # 获取GSE ID（从文件名推断）
        gse_id = os.path.basename(xml_file).replace('.xml', '').replace('_page', '')

        # 存储结果
        result = {
            'GSE_ID': gse_id,
            'XML_File': os.path.basename(xml_file)
        }

        # 对每个问题提问
        for question in questions:
            logger.info(f"  提问: {question}")
            answer = self.ask_question(context, question)
            # 使用简短的关键作为列名
            key = question.lower().replace(' ', '_').replace('?', '').replace('the_', '')
            result[key] = answer
            logger.info(f"  回答: {answer[:100]}...")

        return result

def main():
    # 配置参数
    GEO_DATA_PATH = "/fse/DC/sc_geo_gse"
    MODEL_PATH = "external_inputs/DeepSeek-V3.1"  # 请修改为您的实际模型路径

    # 定义问题
    QUESTIONS = [
        'What is the tissue?',
        'What is the single cell sequencing platform?',
        'What is the disease?'
    ]
    # 查找所有XML文件
    xml_pattern = os.path.join(GEO_DATA_PATH, "**", "*.xml")
    xml_files = glob.glob(xml_pattern, recursive=True)
    if not xml_files:
        logger.error(f"在 {GEO_DATA_PATH} 中未找到XML文件")
        return

    logger.info(f"找到 {len(xml_files)} 个XML文件")
    # 初始化处理器
    processor = DeepSeekProcessor(MODEL_PATH)

    # 存储所有结果
    all_results = []

    # 批量处理文件
    for xml_file in tqdm(xml_files, desc="处理XML文件"):
        try:
            result = processor.process_single_file(xml_file, QUESTIONS)
            if result:
                all_results.append(result)
        except Exception as e:
            logger.error(f"处理文件 {xml_file} 时出错: {e}")
            continue

    # 保存结果到CSV文件
    if all_results:
        output_file = "geo_analysis_results.csv"
        df = pd.DataFrame(all_results)
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        logger.info(f"结果已保存到: {output_file}")

        # 打印摘要
        print("\n=== 处理完成摘要 ===")
        print(f"成功处理文件数: {len(all_results)}")
        print(f"输出文件: {output_file}")

        # 显示前几行结果
        print("\n前3行结果预览:")
        print(df.head(3).to_string(index=False))
    else:
        logger.error("没有成功处理任何文件")

if __name__ == "__main__":
    main()
