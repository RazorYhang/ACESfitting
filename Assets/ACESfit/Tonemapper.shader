Shader "Unlit/Tonemapper"
{
	Properties
	{
		_MainTex ("Texture", 2D) = "white" {}
		_HDR ("HDR multiply", float) = 2
	}
	SubShader
	{
		Tags { "RenderType"="Opaque" }
		LOD 100

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			// make fog work
			#pragma multi_compile_fog
			
			#include "UnityCG.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float2 uv : TEXCOORD0;
				UNITY_FOG_COORDS(1)
				float4 vertex : SV_POSITION;
			};

			sampler2D _MainTex;
			float4 _MainTex_ST;
			float4 param;
			float whitePoint;
			float _HDR;
			
			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = TRANSFORM_TEX(v.uv, _MainTex);
				return o;
			}

			float3 ACESFilmFit_Scalar(float3 x, float4 param, float whitePoint)
			{
				// clamp _x to whitePoint:
				float3 _x = min(x, whitePoint + 0.001);
				// aces curve:
				return (_x*(param.r*_x + param.g)) / (_x*(param.b*_x + param.a) + 1.0f);
			}

			fixed4 frag (v2f i) : SV_Target
			{
				// sample the texture
				fixed4 col = tex2D(_MainTex, i.uv) * _HDR;
				col.rgb = saturate(ACESFilmFit_Scalar(col.rgb, param, whitePoint));
				return col;
			}
			ENDCG
		}
	}
}
