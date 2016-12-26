/* decropper Ver.1.0 */

/* 参考：
 * 「3 スプライン補間」
 * http://www.akita-nct.ac.jp/yamamoto/lecture/2004/5E/interpolation/text/html/node3.html
 * 「1 プログラム方法」
 * http://akita-nct.jp/yamamoto/lecture/2004/5E/interpolation/SplineProgram/html/node1.html
 */

#include <cmath>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "wave.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

bool isLimitNum(const short data, const short lower, const short upper){
	if((data == lower) || (data == upper)) return true;
	return false;
}

short limit(const int a, const short min, const short max){
	if(a < min) return min;
	if(a > max) return max;
	return a;
}

int main(int argc, char *argv[]){
	cout << "[decropper Ver.1.0]" << endl;
	cout << "----------------------------------------" << endl;
	// 標準入力を処理する
	if(argc < 4){
		cout << "usage: decropper wav-file(input) wav-file(output) spline_range" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "input : " << argv[1] << endl;
	cout << "output: " << argv[2] << endl;
	cout << "range : " << argv[3] << endl;
	cout << "----------------------------------------" << endl;
	// ファイルを読み込む
	WAVE wav_file;
	if(!wav_file.load_from_file(argv[1])){
		cout << "error: can't read wav-file!" << endl;
		exit(EXIT_FAILURE);
	}
	WAVE_FORMAT format = wav_file.get_format();
	auto channels = format.num_of_channels;
	auto samples = format.samples_per_sec;
	auto bits = format.bits_per_sample;
	cout << "format id      : " << format.format_id       << endl;
	cout << "num of channels: " << channels << endl;
	cout << "samples per sec: " << samples << endl;
	cout << "bytes per sec  : " << format.bytes_per_sec   << endl;
	cout << "block size     : " << format.block_size      << endl;
	cout << "bits per sample: " << bits << endl;
	cout << "----------------------------------------" << endl;
	vector<vector<short>> data(channels);
	vector<short> min_data(channels), max_data(channels);
	for(auto i = 0; i < channels; ++i){
		// 1チャンネル分を読み込む
		wav_file.get_channel(data[i], i);
		// 最大・最小音を検出する
		min_data[i] = SHRT_MAX;
		max_data[i] = SHRT_MIN;
		for(auto &read_data : data[i]){
			min_data[i] = min_data[i] < read_data ? min_data[i] : read_data;
			max_data[i] = max_data[i] > read_data ? max_data[i] : read_data;
		}
	}
	rsize_t size = data[0].size();
	cout << "data length: " << size << endl;
	cout << "seconds    : " << (static_cast<double>(size) / samples) << endl;
	cout << "----------------------------------------" << endl;
	// 補間処理を行う
	const rsize_t sampling_data_count = abs(std::stoi(argv[3]));
	vector<vector<double>> ip_data(channels, vector<double>(size));
	vector<vector<short>> write_data(channels, vector<short>(size));
	for(auto i = 0; i < channels; ++i){
		cout << "channel" << (i + 1) << ": " << min_data[i] << " ... " << max_data[i] << endl;
		// 補間計算
		for(rsize_t t = 0; t < size; ++t){
			// 最大・最小音でなければそのまま採用する
			if(!isLimitNum(data[i][t], min_data[i], max_data[i])){
				ip_data[i][t] = 1. * data[i][t];
				continue;
			}
			// 最大・最小音に引っかかる場合は、その周辺のデータを収集して判断する
			rsize_t length = 1;
			for(rsize_t t2 = t + 1; t2 < size; ++t2){
				if(isLimitNum(data[i][t2], min_data[i], max_data[i])){
					++length;
				}else{
					break;
				}
			}
			//長さ1なら無視する
			if(length == 1){
				ip_data[i][t] = 1. * data[i][t];
				continue;
			}
			//長さ2以上の場合に補間を行う
			//cout << "[" << (i + 1) << "-" << t << "(" << length << ")]\n";
			{
				// 前後のサンプルを採る
				vector<rsize_t> x;
				vector<short> y;
				rsize_t count = 0;
				for(rsize_t t2 = t - 1; t2 >= 0; --t2){
					if(!isLimitNum(data[i][t2], min_data[i], max_data[i])){
						x.push_back(t2);
						y.push_back(data[i][t2]);
						++count;
						if(count == sampling_data_count) break;
					}
				}
				rsize_t spline_offset = x.size() - 1;
				count = sampling_data_count * 2 - x.size();
				for(rsize_t t2 = t + 1; t2 < size; ++t2){
					if(!isLimitNum(data[i][t2], min_data[i], max_data[i])){
						x.push_back(t2);
						y.push_back(data[i][t2]);
						--count;
						if(count == 0) break;
					}
				}
				// ソートを行う
				for(rsize_t p = 0; p < x.size() - 1; ++p){
					for(rsize_t q = p + 1; q < x.size(); ++q){
						if(x[p] > x[q]){
							rsize_t temp1 = x[p]; x[p] = x[q]; x[q] = temp1;
							short   temp2 = y[p]; y[p] = y[q]; y[q] = temp2;
						}
					}
				}
				//for(rsize_t k = 0; k < x.size(); ++k){
				//	cout << "(" << x[k] << "," << y[k] << "), "; 
				//}
				//cout << "\n";
				// 補間処理を行う
				//計算用行列を構成する
				const rsize_t N = x.size() - 1;
				vector<double> h(N);
				for(rsize_t k = 0; k < N; ++k){
					h[k] = 1. * (x[k + 1] - x[k]);
				}
				vector<double> v(N, 0.0);
				for(rsize_t k = 1; k < N; ++k){
					v[k] = 6. * (y[k + 1] - y[k]) / h[k] - 6. * (y[k] - y[k - 1]) / h[k - 1];
				}
				vector<vector<double>> a(N, vector<double>(N, 0.0));
				for(rsize_t k = 1; k < N; ++k){
					a[k][k] = 2 * (h[k - 1] + h[k]);
				}
				for(rsize_t k = 1; k < N - 1; ++k){
					a[k + 1][k] = h[k];
					a[k][k + 1] = h[k];
				}
				// 連立方程式を解く
				vector<double> u(N + 1, 0.0);
				for(rsize_t k = 1; k < N; ++k){
					// 対角要素で割る
					double div = 1.0 / a[k][k];
					for(rsize_t j = 1; j < N; ++j){
						a[k][j] *= div;
					}
					v[k] *= div;
					// 掃き出し処理を行う
					for(rsize_t l = k + 1; l < N; ++l){
						double mlt = a[l][k];
						for(rsize_t j = k; j < N; ++j){
							a[l][j] -= mlt * a[k][j];
						}
						v[l] -= mlt * v[k];
					}
				}
				for(rsize_t k = N - 1; k >= 1; --k){
					u[k] = v[k];
					for(rsize_t j = k - 1; j >= 1; --j){
						v[j] -= a[j][k] * u[k];
					}
				}
				// 補間を行う
				rsize_t idx = spline_offset;
				double A = (u[idx + 1] - u[idx]) / 6 / (x[idx + 1] - x[idx]);
				double B = u[idx] / 2;
				double C = 1. * (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]) - 1. / 6 * (x[idx + 1] - x[idx]) * (2 * u[idx] + u[idx + 1]);
				double D = 1. * y[idx];
				//cout << A << "," << B << "," << C << "," << D << endl;
				for(rsize_t t2 = t; t2 < t + length; ++t2){
					rsize_t diff = t2 - x[idx];
					ip_data[i][t2] = ((A * diff + B) * diff + C) * diff + D;
					//cout << data[i][t] << "->" << ip_data[i][t2] << endl;
				}
				//
				t += length - 1;
			}
		}
	}
	// リミットを掛ける
	switch(bits){
	case 8:
	{
		double max_ip_data = fabs(ip_data[0][0] - 128);
		for(auto i = 0; i < channels; ++i){
			for(auto &read_data : ip_data[i]){
				max_ip_data = max_ip_data > fabs(read_data - 128) ? max_ip_data : fabs(read_data - 128);
			}
		}
		//cout << max_ip_data << endl;
		for(auto i = 0; i < channels; ++i){
			for(rsize_t t = 0; t < size; ++t){
				write_data[i][t] = limit(static_cast<int>((ip_data[i][t] - 128) / max_ip_data * 127 + 128), 0, 255);
			}
		}
	}
		break;
	case 16:
	{
		double max_ip_data = fabs(ip_data[0][0]);
		for(auto i = 0; i < channels; ++i){
			for(auto &read_data : ip_data[i]){
				max_ip_data = max_ip_data > fabs(read_data) ? max_ip_data : fabs(read_data);
			}
		}
		cout << max_ip_data << endl;
		for(auto i = 0; i < channels; ++i){
			for(rsize_t t = 0; t < size; ++t){
				write_data[i][t] = limit(static_cast<int>(ip_data[i][t] / max_ip_data * 32767), -32768, 32767);
				//if((t >= 9504880) && (t <= 9504890) && (i == 1)){
					//cout << ip_data[i][t] << "->" << write_data[i][t] << "\n";
				//}
			}
		}
	}
		break;
	}
	cout << "----------------------------------------" << endl;
	if(channels == 1){
		wav_file.set_channel(write_data[0], format);
	}else{
		wav_file.set_channel(write_data[0], write_data[1], format);
	}
	wav_file.save_to_file(argv[2]);
}

/*
cd /d F:\ソフトウェア\マルチメディア\音楽\・分析\海苔音源補間計画\第2世代
cl decropper.cpp wave.cpp /O2 /GL /EHsc /nologo
decropper "sine(10Hz)+6dB.wav" "repair.wav" 2
decropper "only my railgun.wav" "repair2.wav" 5
cls
*/
