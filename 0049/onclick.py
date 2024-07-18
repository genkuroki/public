# 以下はChatGPTへの
# 
# * 散布図のデータをマウスのクリックで入力する方法を教えて下さい。
# * 最後のprintの代わりにdataをファイルに保存するにはどうすればよいのですか？
# 
# という質問への回答にあったPythonスクリプトである。
#
# python onclick.py を実行し、
# マウスでクリックしてドットを入力してその窓を閉じると、
# data.csv に入力したデータが保存される。

import matplotlib.pyplot as plt
import csv

data = []

def onclick(event):
    if event.inaxes is not None:
        data.append((event.xdata, event.ydata))
        ax.scatter(event.xdata, event.ydata, c='blue')
        plt.draw()

fig, ax = plt.subplots()
ax.set_title('Click to add points')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()

print(data)

with open('data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y'])
    writer.writerows(data)
