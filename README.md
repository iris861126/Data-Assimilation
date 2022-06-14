# 資料同化(Data Assimilation)

## 目的
使用 Lorenz 40-variable model (Lorenz-96 model)，實作 OI、3DVAR、 EKF 及 EnKF 等資料同化方法，並探討背景誤差協方差矩陣(background error covariance, BEC)及觀測資料分布對於同化結果之影響。

## 模式與資料同化系統之設定
### 模式設定
><img width="1225" alt="image" src="https://user-images.githubusercontent.com/52215353/173612539-911e5a00-03f7-42f0-9d85-6d17ceceb301.png">
>初始化參數設定

><img width="1225" alt="image" src="https://user-images.githubusercontent.com/52215353/173612572-01d3de75-f08d-476c-aa69-f887a0989029.png">
>形成真值、觀測資料及初始分析場之流程

### 資料同化前置處理
><img width="958" alt="image" src="https://user-images.githubusercontent.com/52215353/173625966-dde89405-f366-4602-91b8-bf3c0011c0b9.png">
>穩態 BEC 產生流程圖(Lien, 2020)

><img width="1007" alt="image" src="https://user-images.githubusercontent.com/52215353/173626016-b587200b-176e-4ca5-8ba6-9081f490da0e.png">
>NMC 方法(Lien, 2020)

$$P^b = \alpha E \left(\left[{x_f}^{(t_{2})}-{x_f}^{(t_{1})}\right]\left[{x_f}^{(t_{2})}-{x_f}^{(t_{1})}\right]^T\right) $$

其中，$t_{2}=48hr$、$t_{1}=24hr$、$\alpha$為Rescaling Factor。

NMC法使用多次的預報差值取平均，並加以同質(Homogeneous)及等向 (Isotropic)之假設，得出具有氣候特性的BEC。
><img width="440" alt="image" src="https://user-images.githubusercontent.com/52215353/173618090-1cd2ebaa-762f-4d9e-a36f-258b2112a0dc.png">
>NMC法建立之BEC

>詳細實驗設計及結果部分礙於隱私不公開於此。
