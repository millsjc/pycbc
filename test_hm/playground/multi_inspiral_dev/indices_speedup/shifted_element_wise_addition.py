from manim import *



class ShiftedElementWiseAddition(Scene):
    def add_numbers(self, a, b, result_text_list):
        highlight_box_a = SurroundingRectangle(a, color=YELLOW, buff=0.1)
        highlight_box_b = SurroundingRectangle(b, color=YELLOW, buff=0.1)
        self.play(Create(highlight_box_a), Create(highlight_box_b))

        result = int(a.get_tex_string()) + int(b.get_tex_string())
        result_text = Tex(str(result), color=BLACK)
        result_text.next_to(result_text_list[-1], RIGHT) if result_text_list else result_text.shift(2 * DOWN)
        self.play(Write(result_text))

        self.play(FadeOut(highlight_box_a), FadeOut(highlight_box_b))

        result_text_list.append(result_text)
    
    def construct(self):
        A = [1, 2, 3, 3, 4, 1]
        B = [4, 5, 6, 7]

        A_text = [Tex(str(a), color=BLACK) for a in A]
        B_text = [Tex(str(b), color=BLACK) for b in B]

        for i, a in enumerate(A_text):
            a.shift(RIGHT * i + UP)
        for i, b in enumerate(B_text[::-1]):
            b.shift(LEFT * i)

        self.play(*[Write(a) for a in A_text])
        self.play(*[Write(b) for b in B_text])

        result_text_list = []
        
        # self.add_numbers(a, b, result_text_list)
        
        for i, shift in enumerate(range(-(len(A)-1), len(B)+1)):
            start_idx_A = max(0, shift)
            start_idx_B = max(0, -shift)
            for a, b in zip(A_text[start_idx_A:], B_text[start_idx_B:]):
                self.add_numbers(a, b, result_text_list)

            if i not in [0, 1] and shift != len(B):
                self.play(*[b.animate.shift(RIGHT) for b in B_text])

config.media_width = "60%"
config.background_color = WHITE
config.quality = "high_quality"
