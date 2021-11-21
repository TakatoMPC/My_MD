/*---------------------------------------------------------------------------*/

/*!  \file  MD_main.c
     \brief  calculate forces, positions, and velocities by MD.  */

/*---------------------------------------------------------------------------*/

/* Copyright (C) 2021 Takato ISHIDA  */

/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "MT.h"

#define PI 3.141592f
#define NN 80
#define NRANMX 50000

/*---------------------------------------------------------------------------*/


void Rancal(double *Ran); // 乱数生成関数
void Iniposit(double *rx0, double *ry0, double *Ran, int N0, double L, int NRAN); //初期位置を与える関数
void Inivel(double *velx0, double *vely0, double *Ran, int N0, int NA, double K, double T, int NRAN); //初期速度を与える関数
void Force(double *rx, double *ry, double *fx, double *fy, int N0, double L, double rc, int Switch); //力を計算する関数
void Posit(double *rx, double *ry, double *rx0, double *ry0, double *velx, double *vely, double *fx, double *fy, int N0, int NA, double K, double dt, double L, int Switch); //位置を更新する関数
void Velo(double *velx, double *vely, double *velx0, double *vely0, double *fx, double *fy, double *fx0, double *fy0, int N0, int NA, double K, double dt); //速度を更新する関数
void Store_paramset(FILE *param, int N0, int NA, int NB, double T, double K, double L, double rc, int Switch, double dt);
void Store_state(FILE *res, double *rx, double *ry, double *velx, double *vely, int N0, double time);

int main(void)
    {
          /*  配列の確保
              rx0, ry0: 粒子の初期位置ベクトルの各成分
              rx, ry: 粒子の位置ベクトルの各成分
              velx, vely: 粒子の速度ベクトルの各成分
              fx, fy: 粒子に働く力ベクトルの各成分
              Ran: 一様乱数列が格納してある
          */

    double rx0[NN], ry0[NN], rx[NN], ry[NN], rx00[NN], ry00[NN];
    double velx[NN], vely[NN], velx0[NN], vely0[NN], fx[NN], fy[NN], fx0[NN], fy0[NN];

    double Ran[NRANMX];

          /*  変数宣言
              rxi, ryi: 位置計算のためのバッファ?
              cc0, cc1: バッファ?

              Ntime: 時間積分ループのダミー変数
              Np: 初期データ書き出し用(いらんかも)
              Nopt: 途中のデータ書き出し用(いらんかも)
          */

    double rxi, ryi, Ta, Tb, Ra, Rb;
    double va2 = 0;
    double vb2 = 0;
    int i, ti, count;


    /*Simulation parameters for systems*/
              /*  変数宣言
              T: 無次元化温度
              K: 質量比
              NA: 軽い方A種の数
              NB: 重い方B種の数
              N: NA, NB合計，総数
              dt: 時間刻み
              rc: カットオフ距離
              dens: 粒子数密度
              time: 時間
              L: 系のサイズ(N, 数密度に従属)
              dtsq: dtの二乗，時間積分アルゴリズムで使う
          */
    double T = 1;
    double K = 20;
    int NA = 40;
    int NB = 40;
    // int NA = 20;
    // int NB = 20;
    int N0 = NA + NB;
    double dt = 0.001;
    double rc = 3;
    double dens = 0.1;
    double time = 0;

    double L = sqrt(N0/dens);
    double dtsq = dt*dt;

    /*Simulation parameters for simulation*/
              /*  変数宣言
              tiMAX: 時間ステップの停止までの最大回数
              NRAN: 乱数を配列から取り出す時のインデックス(1~Nをループ)
              Switch: 周期的境界条件 =0 (適用する), otherwise (適用しない)
          */
    int tiMAX = 10000;
    int NRAN = 0;
    int Switch = 0;

    //パラメータセットとinitialと10個の時点の位置・速度をファイル保存する
    FILE *param, *temps, *dif_dis, *res[11];
      param = fopen("Setting.csv", "w");
      temps = fopen("Temperatures.csv", "w");
      dif_dis = fopen("Dif_dist.csv", "w");
      res[0] = fopen("Result0.csv", "w");
      res[1] = fopen("Result1.csv", "w");
      res[2] = fopen("Result2.csv", "w");
      res[3] = fopen("Result3.csv", "w");
      res[4] = fopen("Result4.csv", "w");
      res[5] = fopen("Result5.csv", "w");
      res[6] = fopen("Result6.csv", "w");
      res[7] = fopen("Result7.csv", "w");
      res[8] = fopen("Result8.csv", "w");
      res[9] = fopen("Result9.csv", "w");
      res[10] = fopen("Result10.csv", "w");

    Rancal(Ran);

    //Generate initial configuration
    Iniposit(rx0, ry0, Ran, N0, L, NRAN);
    Inivel(velx0, vely0, Ran, N0, NA, K, T, NRAN);
    Force(rx0, ry0, fx0, fy0, N0, L, rc, Switch);

    //この時点でr0, v0, f0に初期のデータが格納されている．


    for(i=0;i<N0;i++){
      rx[i] = rx0[i];
      ry[i] = ry0[i];
      rx00[i] = rx0[i];
      ry00[i] = ry0[i];
      velx[i] = velx0[i];
      vely[i] = vely0[i];
      fx[i] = fx0[i];
      fy[i] = fy0[i];
      if(i < NA){
          va2 += velx[i]*velx[i] + vely[i]*vely[i];
          }else{
          vb2 += velx[i]*velx[i] + vely[i]*vely[i];
          }
    }
          //分子速度から系の温度を評価してみる
      Ta = 0.5*va2/NA;
      Tb = 0.5*K*vb2/NB;
      //二乗和速度の初期化
      va2 = 0;
      vb2 = 0;
      fprintf(temps, "Ta, Tb\n");
      fprintf(temps, "%f, %f\n", Ta, Tb);
      Ra = 0;
      Rb = 0;
      fprintf(dif_dis, "Ra, Rb\n");
      fprintf(dif_dis, "%f, %f\n", Ra, Rb);




    ////初期・シミュレーション条件をファイル保存する
    Store_paramset(param, N0, NA, NB, T, K, L, rc, Switch, dt);
    Store_state(res[0], rx, ry, velx, vely, N0, time);

    //ここからメインループで時間積分をしていく//
    Switch = 1;
    count = 1;
    for(ti=1;ti<=tiMAX;ti++){
      //速度Verlet法の依存関係を考えて位置，力，速度の順で計算
      Posit(rx, ry, rx0, ry0, velx0, vely0, fx0, fy0, N0, NA, K, dt, L, Switch);
      Force(rx, ry, fx, fy, N0, L, rc, Switch);
      Velo(velx, vely, velx0, vely0, fx, fy, fx0, fy0, N0, NA, K, dt);
      //添字がゼロのものを更新しておく
        for(i=0;i<N0;i++){
          rx0[i] = rx[i];
          ry0[i] = ry[i];
          velx0[i] = velx[i];
          vely0[i] = vely[i];
          fx0[i] = fx[i];
          fy0[i] = fy[i];
          if(i < NA){
            va2 += velx[i]*velx[i] + vely[i]*vely[i];
            Ra += (rx[i]-rx00[i])*(rx[i]-rx00[i])+(rx[i]-rx00[i])*(rx[i]-rx00[i]);
          }else{
            vb2 += velx[i]*velx[i] + vely[i]*vely[i];
            Rb += (rx[i]-rx00[i])*(rx[i]-rx00[i])+(rx[i]-rx00[i])*(rx[i]-rx00[i]);
          }
        }
      //分子速度から系の温度を評価してみる
      Ta = 0.5*va2/NA;
      Tb = 0.5*K*vb2/NB;
      fprintf(temps, "%f, %f\n", Ta, Tb);
      //拡散距離の平均を評価する
      Ra = sqrt(Ra/NA);
      Rb = sqrt(Ra/NB);
      fprintf(dif_dis, "%f, %f\n", Ra, Rb);

      //二乗和速度の初期化
      va2 = 0;
      vb2 = 0;
      //拡散距離の初期化
      Ra = 0;
      Rb = 0;

      if(ti%(tiMAX/10) == 0){
        time = ti*dt;
        Store_state(res[count], rx, ry, velx, vely, N0, time);

        printf("Ta (t=%d) = %f\n", ti, Ta);
        printf("Tb (t=%d) = %f\n", ti, Tb);
        // printf("r1x (t=%d) = %f\t", ti, rx[1]);
        // printf("v1x (t=%d) = %f\n", ti, velx[1]);

        count++;
      }
    }
    fclose(temps);
    fclose(dif_dis);
    return 0;
}

//乱数生成関数Rancal
void Rancal(double *Ran){
  int i;
	init_genrand((unsigned)time(NULL)); //初期化処理
	for(i=0;i<NRANMX;i++){
    Ran[i] = genrand_real1(); //real1のオプションが実数[0,1]の乱数生成処理
	}
}

//初期位置設定関数Iniposit
void Iniposit(double *rx0, double *ry0, double *Ran, int N0,double L, int NRAN){
  int i, j;
  double crx0, cry0, rxij, ryij, rijsq;

  for(i=0;i<N0;i++){
    loop:
    crx0 = L*(Ran[NRAN]-0.5);
    NRAN++;
    cry0 = L*(Ran[NRAN]-0.5);
    NRAN++;

    // printf("Ran = %f\n", Ran[NRAN]);

    if(i!=0){
      for(j=0;j<(i-1);j++){
        rxij = crx0 - rx0[j];
        ryij = cry0 - ry0[j];
        rijsq = rxij*rxij + ryij*ryij;
        //printf("Rsq = %f\n", rijsq);
         if(rijsq<1){
           goto loop; //距離が1以下のところに別のものがある場合はやり直し
        }
      }
    rx0[i] = crx0;
    ry0[i] = cry0;
  }
  }
}

//初期速度設定関数Inivel
void Inivel(double *velx0, double *vely0, double *Ran, int N0, int NA, double K, double T, int NRAN){
  int i;
  double pxa, pya, pxb, pyb;
  pxa = 0;
  pya = 0;
  pxb = 0;
  pyb = 0;

  //マクスウェル分布から一様乱数を使ってサンプリング
  //系の運動量をゼロにするための補正をかける, 4つのpは運動量総和評価のための入れ物

  for(i=0;i<N0;i++){
    if(i < NA){
      velx0[i] = sqrt(-2*T*log(Ran[NRAN]))*cos(2*PI*Ran[NRAN+1]);
      NRAN += 2;
      vely0[i] = sqrt(-2*T*log(Ran[NRAN]))*cos(2*PI*Ran[NRAN+1]);
      NRAN += 2;
      pxa += velx0[i];
      pya += vely0[i];
    }else{
      velx0[i] = sqrt(-2*T/K*log(Ran[NRAN]))*cos(2*PI*Ran[NRAN+1]);
      NRAN += 2;
      vely0[i] = sqrt(-2*T/K*log(Ran[NRAN]))*cos(2*PI*Ran[NRAN+1]);
      NRAN += 2;
      pxb += velx0[i];
      pyb += vely0[i];
    }
	}
    //過剰運動量を1分子あたりに分配
    pxa /= NA;
    pya /= NA;
    pxb /= (N0-NA);
    pyb /= (N0-NA);

    //過剰分を差し引きして速度の調整を行う
    //ここをコメントアウトして調整なしと考察比較してもいいかも
    for(i=0;i<N0;i++){
    if(i < NA){
      velx0[i] -= pxa;
      vely0[i] -= pya;
    }else{
      velx0[i] -= pxb;
      vely0[i] -= pyb;
    }
	}
}

//力計算の関数，主目的: 配列fx[i], fy[i]の更新
void Force(double *rx, double *ry, double *fx, double *fy, int N0, double L, double rc, int Switch){
  int i, j;
  double Rcsq, rxi, ryi, fxi, fyi, rxij, ryij, rijsq; //一時格納のための変数たち
  double fij, fxij, fyij; //力計算のためのバッファ

  Rcsq = rc*rc;

  //Force配列の初期化
  for(i=0;i<N0;i++){
    fx[i] = 0;
    fy[i] = 0;
  }

  for(i=0;i<(N0-1);i++){
    fxi = fx[i];
    fyi = fy[i];

    for(j=(i+1);j<N0;j++){
      //距離の評価
      rxij = rx[i] - rx[j];
      ryij = ry[i] - ry[j];
      // printf("rxij = %f\n", rxij);

      //Switchが0の場合は周期的境界条件の処理を行う
      if(Switch == 0){
        rxij = rxij - rint(rxij/L)*L;
        ryij = ryij - rint(ryij/L)*L;
      }
      //カットオフ外の場合は直ちにループを抜ける
      //計算量削減のために成分レベルでrcより離れてたら計算する必要なし
      if(fabs(rxij) > rc){
        break;
      }
      if(fabs(ryij) > rc){
        break;
      }
      // //ここでノルムでも評価してやる
      rijsq = rxij*rxij + ryij*ryij;
      if(rijsq > Rcsq){
        break;
      }
      //LJp由来の力
      fij = 24*(2*(1/pow(rijsq, 6.0)) - (1/pow(rijsq, 3.0)))/rijsq;
      //printf("rijsq = %f\n", rijsq);
      fxi += (fij*rxij);
      fyi += (fij*ryij);

      //作用反作用を考慮するために相方のjの方の力成分は逆方向のベクトルなので差し引いておく
      fx[j] -= fij*rxij;
      fy[j] -= fij*ryij;
      }

      fx[i] = fxi;
      fy[i] = fyi;




    //確認
    // printf("fx[%d] = %f\n", i, fx[i]);
    // printf("fy0[%d] = %f\n", i, fy[i]);
    }
  }


//位置を更新する関数Posit，目的: 配列rx[i], ry[i]の更新
void Posit(double *rx, double *ry, double *rx0, double *ry0, double *velx, double *vely, double *fx, double *fy, int N0, int NA, double K, double dt, double L, int Switch){
  int i;
  double dt2by2;

  dt2by2 = dt*dt/2;
  for(i=0;i<N0;i++){
    if(i < NA){
      rx[i] = rx0[i] + velx[i]*dt + fx[i]*dt2by2;
      ry[i] = ry0[i] + vely[i]*dt + fy[i]*dt2by2;

  }else{
      rx[i] = rx0[i] + velx[i]*dt + fx[i]*dt2by2/K;
      ry[i] = ry0[i] + vely[i]*dt + fy[i]*dt2by2/K;
  }
  if(Switch == 0){
    rx[i] -= rint(rx[i]/L)*L;
    ry[i] -= rint(ry[i]/L)*L;
  }
  }
}

//速度を更新する関数Velo，目的: 配列velx[i], vely[i]の更新
void Velo(double *velx, double *vely, double *velx0, double *vely0, double *fx, double *fy, double *fx0, double *fy0, int N0, int NA, double K, double dt){
  int i;
  for(i=0;i<N0;i++){
    if(i < NA){
      velx[i] = velx0[i] + 0.5*dt*(fx[i]+fx0[i]);
      vely[i] = vely0[i] + 0.5*dt*(fy[i]+fy0[i]);
  }else{
      velx[i] = velx0[i] + 0.5*dt*(fx[i]+fx0[i])/K;
      vely[i] = vely0[i] + 0.5*dt*(fy[i]+fy0[i])/K;
  }
  }
}




//シミュレーションに用いたパラメータをファイル保存する
void Store_paramset(FILE *param, int N0, int NA, int NB, double T, double K, double L, double rc, int Switch, double dt){
  // fprintf(param, "Used parameter set\n");
  // fprintf(param, "------------------\n");
  fprintf(param, "N, NA, NB, T, K, L, rc, Switch, dt\n");
  fprintf(param, "%d, %d, %d, %f, %f, %f, %f, %d, %f\n", N0, NA, NB, T, K, L, rc, Switch, dt);
  fclose(param);
}

void Store_state(FILE *res, double *rx, double *ry, double *velx, double *vely, int N0, double time){
  int i;
  // fprintf(res, "State of system at t = %f\n", time);
  // fprintf(res, "------------------------\n");
  fprintf(res, "rx, ry, vx, vy\n");

  for(i=0;i<N0;i++){
    fprintf(res, "%f, %f, %f, %f\n", rx[i], ry[i], velx[i], vely[i]);
  }
  fclose(res);
}
