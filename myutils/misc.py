import matplotlib.pyplot as plt


def draw_curve(x_train, x_valid, title="accuary curve", y_label='accuary(%)', save_path=None, legents=[]):
    # 画训练曲线图
    plt.clf()

    x_train, = plt.plot(x_train, linewidth=2, label=legents[0])
    x_valid, = plt.plot(x_valid, linewidth=2, label=legents[1])

    plt.legend()
    plt.title(title, fontsize=12)
    plt.xlabel('epoch', fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    plt.savefig(save_path)


def draw_curve_single(x_train, title="accuary curve", y_label='accuary(%)', save_path=None, legents=[]):
    # 画训练曲线图
    plt.clf()

    x_train, = plt.plot(x_train, linewidth=2, label=legents[0])

    plt.legend()
    plt.title(title, fontsize=12)
    plt.xlabel('epoch', fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    plt.savefig(save_path)


class AverageMeter:
    """ 记录每次训练的损失
    Computes and stores the average and current value
    """
    def __init__(self):
        self.clear()
        self.lst = []
    
    def clear(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0
    
    def clear_all(self):
        self.__init__()
    
    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count
        self.lst.append(self.avg)


def show_cfgs(cfg, trace_func=print):
    for i, k in enumerate(cfg.__dict__):
        trace_func("{} : {}".format(k, cfg.__dict__[k]))