---
title: Git使用方法
tags: [Git]
categories: Git
index_img: /img/default.png
banner_img: /img/default.png
---

git使用方法

<!-- more -->

# Git使用方法

git：分布式版本控制系统

## Git与GIithub远程仓库连接（配置SSH）

（ssh是一种网络协议，用于计算机之间的加密登录）

以下为成功安装git后的操作

### 配置注册用户名和注册邮箱

```
git config --global user.name "xxx"
git config --global user.email "xxx.com"
```

### 查看配置信息

```
git config --list
```

### 生成SSH

```
ssh-keygen -t rsa -C "邮箱"
```

### 获取公钥

SSH文件存放于C:/User/用户/.ssh，其中id_rsa为私钥，id_rsa.pub为公钥。

复制公钥，在GitHub->Setting->SSH and GPG keys->New SSH key，创建一个新的SSH key。

测试SSH连接

```
ssh -T git@github.com
```

## 推送文章到远程仓库

### 在Github上建立新仓库

“+”->New respository

### 建立本地仓库

在目录中新建一个文件夹Git，右键选择Git bash here，执行命令

```
git init
```

将Git文件夹初始化为一个仓库，此时Git文件夹下会出现一个隐藏的.git文件夹。

### 将远程仓库克隆到本地

```
git clone git@github.com:用户名/仓库名.git
```

### 管理仓库中文件

### 命令add

在仓库文件夹下git bash here，执行命令

```
git add 单个文件
git add 文件夹1/ 文件夹2/ ...
git add .
```

### 命令commit

```
git commit -m "注释"
```

### 命令push

```
git push -u origin main
```

至此，便完成了“远程仓库的建立->本地仓库克隆->管理仓库->推送本地仓库到远程”的过程。



