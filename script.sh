#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i u0.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0256_25x.mp4
#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrl.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrl256_25x.mp4
#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i u0.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0256_128y.mp4
#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrl.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrl256_128y.mp4
#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i u0.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0256_247z.mp4
#ffmpeg -r 4 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrl.%03d.pgm -vframes 500 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrl256_247z.mp4
#ffmpeg -r 4 -f image2 -s 1024x1024 -timecode_frame_start 0 -i ctestx.%04d.pgm -vframes 501 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrl512_239x.mp4
ffmpeg -r 4 -f image2 -s 1024x1024 -timecode_frame_start 0 -i u0.%04d.pgm -vframes 151 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0512_256y_2.mp4
#ffmpeg -r 2 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrlx.%03d.pgm -vframes 401 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrlx32.mp4
#ffmpeg -r 2 -f image2 -s 512x512 -timecode_frame_start 0 -i u0y.%03d.pgm -vframes 101 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0y128.mp4
#ffmpeg -r 2 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrly.%03d.pgm -vframes 401 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrly128.mp4
#ffmpeg -r 2 -f image2 -s 512x512 -timecode_frame_start 0 -i u0z.%03d.pgm -vframes 101 -vcodec libx264 -crf 25  -pix_fmt yuv420p u0z246.mp4
#ffmpeg -r 2 -f image2 -s 512x512 -timecode_frame_start 0 -i ucrlz.%03d.pgm -vframes 401 -vcodec libx264 -crf 25  -pix_fmt yuv420p ucrlz242.mp4
#for QIAOYIN in {100..250}
#do
#   QIAOYIN2=$QIAOYIN'.pgm'
#   rm -f utestx.$QIAOYIN2
#   rm -f ctest.$QIAOYIN2
#   rm -f utesty.$QIAOYIN2
#   rm -f ctesty.$QIAOYIN2
#   rm -f utestz.$QIAOYIN2
#   rm -f ctestz.$QIAOYIN2
#done
#echo 'done'
